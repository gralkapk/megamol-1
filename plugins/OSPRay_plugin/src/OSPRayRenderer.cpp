/*
 * OSPRayRenderer.cpp
 * Copyright (C) 2009-2017 by MegaMol Team
 * Alle Rechte vorbehalten.
 */

#include "stdafx.h"
#include "OSPRayRenderer.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#include "vislib/graphics/CameraParamsStore.h"
#include "vislib/graphics/gl/IncludeAllGL.h"
#include "vislib/graphics/gl/ShaderSource.h"
#include "vislib/math/Vector.h"
#include "mmcore/utility/log/Log.h"

#include "mmcore/CoreInstance.h"

#include <chrono>
#include <functional>

#include "ospray/ospray_cpp.h"

#include <sstream>
#include <stdint.h>
#include <corecrt_math_defines.h>

using namespace megamol::ospray;

/*
ospray::OSPRayRenderer::OSPRaySphereRenderer
*/
OSPRayRenderer::OSPRayRenderer(void)
    : AbstractOSPRayRenderer()
	, _cam()
    , _osprayShader()
    , _getStructureSlot("getStructure", "Connects to an OSPRay structure")

{
    this->_getStructureSlot.SetCompatibleCall<CallOSPRayStructureDescription>();
    this->MakeSlotAvailable(&this->_getStructureSlot);

    _imgSize = {0,0};
    _time = 0;
    _framebuffer = nullptr;
    _renderer = nullptr;
    _camera = nullptr;
    _world = nullptr;

    _accum_time.count = 0;
    _accum_time.amount = 0;
}


/*
ospray::OSPRayRenderer::~OSPRaySphereRenderer
*/
OSPRayRenderer::~OSPRayRenderer(void) {
    this->_osprayShader.Release();
    this->Release();
}


/*
ospray::OSPRayRenderer::create
*/
bool OSPRayRenderer::create() {
    ASSERT(IsAvailable());

    vislib::graphics::gl::ShaderSource vert, frag;

    if (!instance()->ShaderSourceFactory().MakeShaderSource("ospray::vertex", vert)) {
        return false;
    }
    if (!instance()->ShaderSourceFactory().MakeShaderSource("ospray::fragment", frag)) {
        return false;
    }

    try {
        if (!this->_osprayShader.Create(vert.Code(), vert.Count(), frag.Code(), frag.Count())) {
            megamol::core::utility::log::Log::DefaultLog.WriteMsg(
                megamol::core::utility::log::Log::LEVEL_ERROR, "Unable to compile ospray shader: Unknown error\n");
            return false;
        }
    } catch (vislib::graphics::gl::AbstractOpenGLShader::CompileException ce) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(megamol::core::utility::log::Log::LEVEL_ERROR,
            "Unable to compile ospray shader: (@%s): %s\n",
            vislib::graphics::gl::AbstractOpenGLShader::CompileException::CompileActionName(ce.FailedAction()),
            ce.GetMsgA());
        return false;
    } catch (vislib::Exception e) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(
            megamol::core::utility::log::Log::LEVEL_ERROR, "Unable to compile ospray shader: %s\n", e.GetMsgA());
        return false;
    } catch (...) {
        megamol::core::utility::log::Log::DefaultLog.WriteMsg(
            megamol::core::utility::log::Log::LEVEL_ERROR, "Unable to compile ospray shader: Unknown exception\n");
        return false;
    }

    // this->initOSPRay(device);
    this->setupTextureScreen();
    // this->setupOSPRay(renderer, camera, world, "scivis");#

    return true;
}

/*
ospray::OSPRayRenderer::release
*/
void OSPRayRenderer::release() {

    ospShutdown();

    releaseTextureScreen();

}

/*
ospray::OSPRayRenderer::Render
*/
bool OSPRayRenderer::Render(megamol::core::view::CallRender3D_2& cr) {
    this->initOSPRay();

    
    // if user wants to switch renderer
    if (this->_rd_type.IsDirty()) {
        //ospRelease(_camera);
        //ospRelease(_world);
        //ospRelease(_renderer);
        switch (this->_rd_type.Param<core::param::EnumParam>()->Value()) {
        case PATHTRACER:
            this->setupOSPRay("pathtracer");
            this->_rd_type_string = "pathtracer";
            break;
        case MPI_RAYCAST: //< TODO: Probably only valid if device is a "mpi_distributed" device
            this->setupOSPRay("mpi_raycast");
            this->_rd_type_string = "mpi_raycast";
            break;
        default:
            this->setupOSPRay("scivis");
            this->_rd_type_string = "scivis";
        }
        _renderer_has_changed = true;
        this->_rd_type.ResetDirty();
    }

    if (&cr == nullptr) return false;


    CallOSPRayStructure* os = this->_getStructureSlot.CallAs<CallOSPRayStructure>();
    if (os == nullptr) return false;
    // read data
    os->setStructureMap(&_structureMap);
    os->setTime(cr.Time());
    if (!os->fillStructureMap()) return false;
    // check if data has changed
    _data_has_changed = false;
    _material_has_changed = false;
    _transformation_has_changed = false;
    for (auto element : this->_structureMap) {
        auto structure = element.second;
        if (structure.dataChanged) {
            _data_has_changed = true;
        }
        if (structure.materialChanged) {
            _material_has_changed = true;
        }
        if (structure.transformationChanged) {
            _transformation_has_changed = true;
        }
    }

    // Light setup
    _light_has_changed = this->GetLights();

    core::view::Camera_2 tmp_newcam;
    cr.GetCamera(tmp_newcam);
    cam_type::snapshot_type snapshot;
    cam_type::matrix_type viewTemp, projTemp;

	// Generate complete snapshot and calculate matrices
    tmp_newcam.calc_matrices(snapshot, viewTemp, projTemp, core::thecam::snapshot_content::all);

    // check data and camera hash
    if (_cam.eye_position().x() != tmp_newcam.eye_position().x() ||
	    _cam.eye_position().y() != tmp_newcam.eye_position().y() ||
	    _cam.eye_position().z() != tmp_newcam.eye_position().z() ||
	    _cam.view_vector() != tmp_newcam.view_vector()
	    ) {
	    _cam_has_changed = true;
    } else {
	    _cam_has_changed = false;
    }

    // Generate complete snapshot and calculate matrices
    _cam = tmp_newcam;

    // glDisable(GL_CULL_FACE);

    // new framebuffer at resize action
    // bool triggered = false;
    if (_imgSize[0] != _cam.resolution_gate().width() ||
        _imgSize[1] != _cam.resolution_gate().height() || _accumulateSlot.IsDirty()) {
        // triggered = true;
        // Breakpoint for Screenshooter debugging
        // if (framebuffer != NULL) ospFreeFrameBuffer(framebuffer);
        _imgSize[0] = _cam.resolution_gate().width();
        _imgSize[1] = _cam.resolution_gate().height();
        _framebuffer = std::make_shared<::ospray::cpp::FrameBuffer>(_imgSize[0], _imgSize[1], OSP_FB_RGBA8, OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
        _db.resize(_imgSize[0] * _imgSize[1]);
        _framebuffer->commit();
    }

    //// if user wants to switch renderer
    // if (this->rd_type.IsDirty()) {
    //    ospRelease(camera);
    //    ospRelease(world);
    //    ospRelease(renderer);
    //    switch (this->rd_type.Param<core::param::EnumParam>()->Value()) {
    //    case PATHTRACER:
    //        this->setupOSPRay(renderer, camera, world, "pathtracer");
    //        break;
    //    case MPI_RAYCAST: //< TODO: Probably only valid if device is a "mpi_distributed" device
    //        this->setupOSPRay(renderer, camera, world, "mpi_raycast");
    //        break;
    //    default:
    //        this->setupOSPRay(renderer, camera, world, "scivis");
    //    }
    //    renderer_has_changed = true;
    //}
    setupOSPRayCamera(_cam);
    _camera->commit();

    _osprayShader.Enable();
    // if nothing changes, the image is rendered multiple times
    if (_data_has_changed || _material_has_changed || _light_has_changed || _cam_has_changed || _renderer_has_changed ||
        _transformation_has_changed || !(this->_accumulateSlot.Param<core::param::BoolParam>()->Value()) ||
        _frameID != static_cast<size_t>(cr.Time()) || this->InterfaceIsDirty()) {

        std::array<float, 4> eyeDir = {
            _cam.view_vector().x(), _cam.view_vector().y(), _cam.view_vector().z(), _cam.view_vector().w()};
        if (_data_has_changed || _frameID != static_cast<size_t>(cr.Time()) || _renderer_has_changed) {
            // || this->InterfaceIsDirty()) {
            if (!this->generateRepresentations()) return false;
            this->createInstances();
            std::vector<::ospray::cpp::Instance> instanceArray;
            std::transform(_instances.begin(), _instances.end(), std::back_inserter(instanceArray), second(_instances));
            _world->setParam("instance", ::ospray::cpp::CopiedData(instanceArray));

            // Enable Lights
            this->fillLightArray(eyeDir);
            _world->setParam("light", ::ospray::cpp::CopiedData(_lightArray));

            // Commiting world and measuring time
            auto t1 = std::chrono::high_resolution_clock::now();
            _world->commit();
            auto t2 = std::chrono::high_resolution_clock::now();
            const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1).count();
            megamol::core::utility::log::Log::DefaultLog.WriteMsg(
                242, "[OSPRayRenderer] Commiting World took: %d microseconds", duration);
        }
        if (_material_has_changed && !_data_has_changed) {
            this->changeMaterial();
        }
        if (_transformation_has_changed  || _material_has_changed && !_data_has_changed) {
            this->changeTransformation();
            std::vector<::ospray::cpp::Instance> instanceArray;
            std::transform(_instances.begin(), _instances.end(), std::back_inserter(instanceArray), second(_instances));
            _world->setParam("instance", ::ospray::cpp::CopiedData(instanceArray));
            _world->commit();
        }
        if (_light_has_changed && !_data_has_changed) {
            this->fillLightArray(eyeDir);
            _world->setParam("light", ::ospray::cpp::CopiedData(_lightArray));
            _world->commit();
        }


        this->InterfaceResetDirty();
        _time = cr.Time();
        _frameID = static_cast<size_t>(cr.Time());
        _renderer_has_changed = false;

        /*
            if (this->maxDepthTexture) {
                ospRelease(this->maxDepthTexture);
            }
            this->maxDepthTexture = getOSPDepthTextureFromOpenGLPerspective(*cr);
        */
        RendererSettings(cr.BackgroundColor());


        if (this->_useDB.Param<core::param::BoolParam>()->Value()) {
            // far distance
            float far_clip = _cam.far_clipping_plane();
            std::vector<float> far_dist(_imgSize[0] * _imgSize[1], far_clip);
            rkcommon::math::vec2i imgSize = {
                _imgSize[0],
                _imgSize[1]
            };

            auto depth_texture_data = ::ospray::cpp::CopiedData(far_dist.data(), OSP_FLOAT, imgSize);
            depth_texture_data.commit();
            auto depth_texture = ::ospray::cpp::Texture("texture2d");
            depth_texture.setParam("format", OSP_TEXTURE_R32F);
            depth_texture.setParam("filter", OSP_TEXTURE_FILTER_NEAREST);
            depth_texture.setParam("data", depth_texture_data);
            depth_texture.commit();

            _renderer->setParam("map_maxDepth", depth_texture);
        }
        _renderer->commit();

        // setup framebuffer and measure time
        auto t1 = std::chrono::high_resolution_clock::now();

        _framebuffer->clear();//(OSP_FB_COLOR | OSP_FB_DEPTH | OSP_FB_ACCUM);
        _framebuffer->renderFrame(*_renderer, *_camera, *_world);

        // get the texture from the framebuffer
        _fb = reinterpret_cast<uint32_t*>(_framebuffer->map(OSP_FB_COLOR));

        auto t2 = std::chrono::high_resolution_clock::now();
        const auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);

        _accum_time.amount += duration.count();
        _accum_time.count += 1;
        if (_accum_time.amount >= static_cast<unsigned long long int>(1e6)) {
            const unsigned long long int mean_rendertime = _accum_time.amount / _accum_time.count;
            megamol::core::utility::log::Log::DefaultLog.WriteMsg(242, "[OSPRayRenderer] Rendering took: %d microseconds", mean_rendertime);
            _accum_time.count = 0;
            _accum_time.amount = 0;
        }


        float* db;
        if (this->_useDB.Param<core::param::BoolParam>()->Value()) {
             db = static_cast<float*>(_framebuffer->map(OSP_FB_DEPTH));
            _db = std::vector<float>(db, db + _imgSize[0] * _imgSize[1]);
        
            getOpenGLDepthFromOSPPerspective(_db.data());
        }

        // write a sequence of single pictures while the screenshooter is running
        // only for debugging
        // if (triggered) {
        //    std::ostringstream oss;
        //    oss << "ospframe" << this->number << ".ppm";
        //    std::string bla = oss.str();
        //    const char* fname = bla.c_str();
        //    osp::vec2i isize;
        //    isize.x = cr->GetCameraParameters()->TileRect().GetSize().GetWidth();
        //    isize.y = cr->GetCameraParameters()->TileRect().GetSize().GetHeight();
        //    writePPM(fname, isize, fb);
        //    this->number++;
        //}

        //std::string fname("blub.ppm");
        //writePPM(fname, _imgSize, _fb);
        

        this->renderTexture2D(_osprayShader, _fb, _db.data(), _imgSize[0], _imgSize[1], cr);

        // clear stuff
         _framebuffer->unmap(_fb);
        if (this->_useDB.Param<core::param::BoolParam>()->Value()) {
            _framebuffer->unmap(db);
        }

        //auto dvce_ = ospGetCurrentDevice();
        //auto error_ = std::string(ospDeviceGetLastErrorMsg(dvce_));
        //megamol::core::utility::log::Log::DefaultLog.WriteError(std::string("OSPRAY last ERROR: " + error_).c_str());

        this->releaseOSPRayStuff();


    } else {
        _framebuffer->renderFrame(*_renderer, *_camera, *_world);
        _fb = reinterpret_cast<uint32_t*>(_framebuffer->map(OSP_FB_COLOR));

        this->renderTexture2D(_osprayShader, _fb, _db.data(), _imgSize[0], _imgSize[1], cr);
        _framebuffer->unmap(_fb);
    }

    _osprayShader.Disable();

    return true;
}

/*
ospray::OSPRayRenderer::InterfaceIsDirty()
*/
bool OSPRayRenderer::InterfaceIsDirty() {
    if (this->AbstractIsDirty()) {
        return true;
    } else {
        return false;
    }
}

/*
ospray::OSPRayRenderer::InterfaceResetDirty()
*/
void OSPRayRenderer::InterfaceResetDirty() { this->AbstractResetDirty(); }


/*
 * ospray::OSPRayRenderer::GetExtents
 */
bool OSPRayRenderer::GetExtents(megamol::core::view::CallRender3D_2& cr) {

    if (&cr == NULL) return false;
    CallOSPRayStructure* os = this->_getStructureSlot.CallAs<CallOSPRayStructure>();
    if (os == NULL) return false;
    os->setTime(static_cast<int>(cr.Time()));
    os->setExtendMap(&(this->_extendMap));
    if (!os->fillExtendMap()) return false;

    megamol::core::BoundingBoxes_2 finalBox;
    unsigned int frameCnt = 0;
    for (auto pair : this->_extendMap) {
        auto element = pair.second;

        if (frameCnt == 0) {
            if (element.boundingBox->IsBoundingBoxValid()) {
                finalBox.SetBoundingBox(element.boundingBox->BoundingBox());
            } else if (element.boundingBox->IsClipBoxValid()) {
                finalBox.SetBoundingBox(element.boundingBox->ClipBox());
            } else {
                finalBox.SetBoundingBox(vislib::math::Cuboid<float>(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f));
            }
            if (element.boundingBox->IsClipBoxValid()) {
                finalBox.SetClipBox(element.boundingBox->ClipBox());
            } else if (element.boundingBox->IsBoundingBoxValid()) {
                finalBox.SetClipBox(element.boundingBox->BoundingBox());
            } else {
                finalBox.SetClipBox(vislib::math::Cuboid<float>(-1.0f, -1.0f, -1.0f, 1.0f, 1.0f, 1.0f));
            }

        } else {
            if (element.boundingBox->IsBoundingBoxValid()) {
                vislib::math::Cuboid<float> box(finalBox.BoundingBox());
                box.Union(element.boundingBox->BoundingBox());
                finalBox.SetBoundingBox(box);
            } else if (element.boundingBox->IsClipBoxValid()) {
                vislib::math::Cuboid<float> box(finalBox.BoundingBox());
                box.Union(element.boundingBox->BoundingBox());
                finalBox.SetBoundingBox(box);
            }
            if (element.boundingBox->IsClipBoxValid()) {
                vislib::math::Cuboid<float> box(finalBox.ClipBox());
                box.Union(element.boundingBox->ClipBox());
                finalBox.SetClipBox(box);
            } else if (element.boundingBox->IsBoundingBoxValid()) {
                vislib::math::Cuboid<float> box(finalBox.ClipBox());
                box.Union(element.boundingBox->BoundingBox());
                finalBox.SetClipBox(box);
            }
        }
        frameCnt = vislib::math::Max(frameCnt, element.timeFramesCount);
    }
    cr.SetTimeFramesCount(frameCnt);

    cr.AccessBoundingBoxes().SetBoundingBox(finalBox.BoundingBox());
    cr.AccessBoundingBoxes().SetBoundingBox(finalBox.ClipBox());

    return true;
}

void OSPRayRenderer::getOpenGLDepthFromOSPPerspective(float* db) {

    const float fovy = _cam.aperture_angle();
    const float aspect = _cam.resolution_gate_aspect();
    const float zNear = _cam.near_clipping_plane();
    const float zFar = _cam.far_clipping_plane();

    const glm::vec3 cameraUp = {_cam.up_vector().x(), _cam.up_vector().y(), _cam.up_vector().z()};
    const glm::vec3 cameraDir = {_cam.view_vector().x(), _cam.view_vector().y(), _cam.view_vector().z()};

    // map OSPRay depth buffer from provided frame buffer
    auto ospDepthBuffer = static_cast<float*>(_framebuffer->map(OSP_FB_DEPTH));

    const auto ospDepthBufferWidth = static_cast<const size_t>(_imgSize[0]);
    const auto ospDepthBufferHeight = static_cast<const size_t>(_imgSize[1]);

    // transform from ray distance t to orthogonal Z depth
    auto dir_du = glm::normalize(glm::cross(cameraDir, cameraUp));
    auto dir_dv = glm::normalize(glm::cross(dir_du, cameraDir));

    const float imagePlaneSizeY = 2.f * tanf(fovy / 2.f * M_PI / 180.f);
    const float imagePlaneSizeX = imagePlaneSizeY * aspect;

    dir_du *= imagePlaneSizeX;
    dir_dv *= imagePlaneSizeY;

    const auto dir_00 = cameraDir - .5f * dir_du - .5f * dir_dv;

    const float A = -(zFar + zNear) / (zFar - zNear);
    const float B = -2. * zFar * zNear / (zFar - zNear);

    int j, i;
#pragma omp parallel for private(i)
    for (j = 0; j < ospDepthBufferHeight; j++) {
        for (i = 0; i < ospDepthBufferWidth; i++) {
            const auto dir_ij = glm::normalize(dir_00 + float(i) / float(ospDepthBufferWidth - 1) * dir_du +
                                                      float(j) / float(ospDepthBufferHeight - 1) * dir_dv);

            const float tmp = ospDepthBuffer[j * ospDepthBufferWidth + i];// * dot(cameraDir, dir_ij);
            float res = 0.5 * (-A * tmp + B) / tmp + 0.5;
            if (!std::isfinite(res)) res = 1.0f;
            db[j * ospDepthBufferWidth + i] = res;
        }
    }
    // unmap OSPRay depth buffer
    _framebuffer->unmap(ospDepthBuffer);
}
