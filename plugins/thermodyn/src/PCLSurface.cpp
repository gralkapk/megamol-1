#include "PCLSurface.h"

#include "geometry_calls/MultiParticleDataCall.h"
#include "mesh/MeshCalls.h"

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl/surface/gp3.h>
#include <pcl/io/obj_io.h>
#include <pcl/surface/mls.h>
#include <pcl/surface/marching_cubes_rbf.h>
#include <pcl/surface/marching_cubes_hoppe.h>
#include <pcl/geometry/triangle_mesh.h>


megamol::thermodyn::PCLSurface::PCLSurface() : out_data_slot_("outData", ""), in_data_slot_("inData", "") {
    out_data_slot_.SetCallback(mesh::CallMesh::ClassName(), mesh::CallMesh::FunctionName(0), &PCLSurface::get_data_cb);
    out_data_slot_.SetCallback(
        mesh::CallMesh::ClassName(), mesh::CallMesh::FunctionName(1), &PCLSurface::get_extent_cb);
    MakeSlotAvailable(&out_data_slot_);

    in_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_data_slot_);
}


megamol::thermodyn::PCLSurface::~PCLSurface() {
    this->Release();
}


bool megamol::thermodyn::PCLSurface::create() {
    return true;
}


void megamol::thermodyn::PCLSurface::release() {}


bool megamol::thermodyn::PCLSurface::get_data_cb(core::Call& c) {
    auto out_call = dynamic_cast<mesh::CallMesh*>(&c);
    if (out_call == nullptr)
        return false;

    auto in_call = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_call == nullptr)
        return false;

    if (!(*in_call)(0))
        return false;

    if (in_call->DataHash() != in_data_hash_ || in_call->FrameID() != frame_id_) {
        auto const pl_count = in_call->GetParticleListCount();

        for (std::remove_const_t<decltype(pl_count)> pl_idx = 0; pl_idx < pl_count; ++pl_idx) {
            auto const& particles = in_call->AccessParticles(pl_idx);

            pcl::PointCloud<pcl::PointXYZ>::Ptr p_cloud(new pcl::PointCloud<pcl::PointXYZ>);
            pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);

            auto const p_count = particles.GetCount();
            auto const x_acc = particles.GetParticleStore().GetXAcc();
            auto const y_acc = particles.GetParticleStore().GetYAcc();
            auto const z_acc = particles.GetParticleStore().GetZAcc();
            auto const dx_acc = particles.GetParticleStore().GetDXAcc();
            auto const dy_acc = particles.GetParticleStore().GetDYAcc();
            auto const dz_acc = particles.GetParticleStore().GetDZAcc();
            for (std::remove_const_t<decltype(p_count)> p_idx = 0; p_idx < p_count; ++p_idx) {
                p_cloud->emplace_back(x_acc->Get_f(p_idx), y_acc->Get_f(p_idx), z_acc->Get_f(p_idx));
                normals->emplace_back(dx_acc->Get_f(p_idx), dy_acc->Get_f(p_idx), dz_acc->Get_f(p_idx));
            }

            /*pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> n;
            pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
            pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
            tree->setInputCloud(p_cloud);
            n.setInputCloud(p_cloud);
            n.setSearchMethod(tree);
            n.setKSearch(20);
            n.compute(*normals);*/

            pcl::PointCloud<pcl::PointNormal>::Ptr cloud_with_normals(new pcl::PointCloud<pcl::PointNormal>);
            pcl::concatenateFields(*p_cloud, *normals, *cloud_with_normals);

            pcl::search::KdTree<pcl::PointNormal>::Ptr tree2(new pcl::search::KdTree<pcl::PointNormal>);
            tree2->setInputCloud(cloud_with_normals);

            pcl::PolygonMesh triangles;

            //pcl::PointCloud<pcl::PointNormal>::Ptr cloud_mls(new pcl::PointCloud<pcl::PointNormal>);
            //pcl::MovingLeastSquares<pcl::PointNormal, pcl::PointNormal> mls;
            ////mls.setComputeNormals(true);
            //// Set parameters
            //mls.setInputCloud(cloud_with_normals);
            //mls.setPolynomialOrder(2);
            //mls.setSearchMethod(tree2);
            //mls.setSearchRadius(6);
            //mls.process(*cloud_mls);

            //pcl::search::KdTree<pcl::PointNormal>::Ptr tree3(new pcl::search::KdTree<pcl::PointNormal>);
            //tree3->setInputCloud(cloud_mls);

            pcl::MarchingCubesHoppe<pcl::PointNormal> mc_rbf;
            mc_rbf.setInputCloud(cloud_with_normals);
            mc_rbf.setSearchMethod(tree2);
            mc_rbf.setIsoLevel(0.5);
            mc_rbf.setGridResolution(128, 128, 128);
            mc_rbf.reconstruct(triangles);


            //pcl::PointCloud<pcl::PointNormal>::Ptr cloud_mls(new pcl::PointCloud<pcl::PointNormal>);
            //pcl::MovingLeastSquares<pcl::PointNormal, pcl::PointNormal> mls;
            ////mls.setComputeNormals(true);
            //// Set parameters
            //mls.setInputCloud(cloud_with_normals);
            //mls.setPolynomialOrder(2);
            //mls.setSearchMethod(tree2);
            //mls.setSearchRadius(6);
            //mls.process(*cloud_mls);

            //pcl::search::KdTree<pcl::PointNormal>::Ptr tree3(new pcl::search::KdTree<pcl::PointNormal>);
            //tree2->setInputCloud(cloud_mls);


            //pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;

            ////gp3.setSearchRadius(0.025);
            //gp3.setSearchRadius(15);

            //gp3.setMu(2.5);
            //gp3.setMaximumNearestNeighbors(512);
            //gp3.setMaximumSurfaceAngle(M_PI / 2); // 45 degrees
            //gp3.setMinimumAngle(M_PI / 18);       // 10 degrees
            //gp3.setMaximumAngle(2 * M_PI / 3);    // 120 degrees
            //gp3.setNormalConsistency(false);

            //gp3.setInputCloud(cloud_mls);
            //gp3.setSearchMethod(tree3);
            //gp3.reconstruct(triangles);

            //std::vector<int> parts = gp3.getPartIDs();
            //std::vector<int> states = gp3.getPointStates();

            //pcl::io::saveOBJFile("test.obj", triangles);


            // convert to mesh
            auto const point_count = triangles.cloud.width * triangles.cloud.height;
            auto const point_size = triangles.cloud.data.size() / point_count;
            auto const face_count = triangles.polygons.size();

            size_t field_map[6] = {0};
            for (decltype(triangles.cloud.fields.size()) f_idx = 0; f_idx < triangles.cloud.fields.size(); ++f_idx) {
                if (triangles.cloud.fields[f_idx].name == "x") {
                    field_map[0] = triangles.cloud.fields[f_idx].offset;
                } else if (triangles.cloud.fields[f_idx].name == "y") {
                    field_map[1] = triangles.cloud.fields[f_idx].offset;
                } else if (triangles.cloud.fields[f_idx].name == "z") {
                    field_map[2] = triangles.cloud.fields[f_idx].offset;
                } else if (triangles.cloud.fields[f_idx].name == "normal_x") {
                    field_map[3] = triangles.cloud.fields[f_idx].offset;
                } else if (triangles.cloud.fields[f_idx].name == "normal_y") {
                    field_map[4] = triangles.cloud.fields[f_idx].offset;
                } else if (triangles.cloud.fields[f_idx].name == "normal_z") {
                    field_map[5] = triangles.cloud.fields[f_idx].offset;
                }
            }

            vertices_.clear();
            vertices_.resize(point_count);
            normals_.clear();
            normals_.resize(point_count);
            for (std::remove_const_t<decltype(point_count)> p_idx = 0; p_idx < point_count; ++p_idx) {
                auto const x = *reinterpret_cast<float*>(&triangles.cloud.data[p_idx * point_size + field_map[0]]);
                auto const y = *reinterpret_cast<float*>(&triangles.cloud.data[p_idx * point_size + field_map[1]]);
                auto const z = *reinterpret_cast<float*>(&triangles.cloud.data[p_idx * point_size + field_map[2]]);
                /*auto const nx = *reinterpret_cast<float*>(&triangles.cloud.data[p_idx * point_size + field_map[3]]);
                auto const ny = *reinterpret_cast<float*>(&triangles.cloud.data[p_idx * point_size + field_map[4]]);
                auto const nz = *reinterpret_cast<float*>(&triangles.cloud.data[p_idx * point_size + field_map[5]]);*/
                vertices_[p_idx] = glm::vec3(x, y, z);
                //normals_[p_idx] = glm::vec3(nx, ny, nz);
            }

            indices_.clear();
            indices_.resize(face_count);
            std::vector<glm::vec3> face_normals(face_count);
            std::vector<std::pair<glm::vec3, unsigned int>> vert_norm_acc(
                point_count, std::make_pair<glm::vec3, unsigned int>(glm::vec3(0), 0));
            for (std::remove_const_t<decltype(face_count)> f_idx = 0; f_idx < face_count; ++f_idx) {
                auto const i0 = triangles.polygons[f_idx].vertices[0];
                auto const i1 = triangles.polygons[f_idx].vertices[1];
                auto const i2 = triangles.polygons[f_idx].vertices[2];
                indices_[f_idx] = glm::uvec3(i0, i1, i2);
                face_normals[f_idx] = glm::normalize(glm::cross((vertices_[i1] - vertices_[i0]), (vertices_[i2] - vertices_[i0])));
                vert_norm_acc[i0].first += face_normals[f_idx];
                vert_norm_acc[i0].second += 1;
                vert_norm_acc[i1].first += face_normals[f_idx];
                vert_norm_acc[i1].second += 1;
                vert_norm_acc[i2].first += face_normals[f_idx];
                vert_norm_acc[i2].second += 1;
            }

            for (std::remove_const_t<decltype(point_count)> p_idx = 0; p_idx < point_count; ++p_idx) {
                auto const norm = vert_norm_acc[p_idx].first / static_cast<float>(vert_norm_acc[p_idx].second);
                normals_[p_idx] = -norm;
            }

            mesh::MeshDataAccessCollection::IndexData idx_data;
            idx_data.data = reinterpret_cast<uint8_t*>(indices_.data());
            idx_data.byte_size = indices_.size() * sizeof(decltype(indices_)::value_type);
            idx_data.type = mesh::MeshDataAccessCollection::UNSIGNED_INT;

            mesh::MeshDataAccessCollection::VertexAttribute vert_attr;
            vert_attr.data = reinterpret_cast<uint8_t*>(vertices_.data());
            vert_attr.component_cnt = 3;
            vert_attr.component_type = mesh::MeshDataAccessCollection::FLOAT;
            vert_attr.offset = 0;
            vert_attr.stride = sizeof(decltype(vertices_)::value_type);
            vert_attr.byte_size = vertices_.size() * sizeof(decltype(vertices_)::value_type);
            vert_attr.semantic = mesh::MeshDataAccessCollection::POSITION;

            mesh::MeshDataAccessCollection::VertexAttribute normal_attr;
            normal_attr.data = reinterpret_cast<uint8_t*>(normals_.data());
            normal_attr.component_cnt = 3;
            normal_attr.component_type = mesh::MeshDataAccessCollection::FLOAT;
            normal_attr.offset = 0;
            normal_attr.stride = sizeof(decltype(normals_)::value_type);
            normal_attr.byte_size = normals_.size() * sizeof(decltype(normals_)::value_type);
            normal_attr.semantic = mesh::MeshDataAccessCollection::NORMAL;

            mesh::MeshDataAccessCollection::Mesh mesh;
            mesh.primitive_type = mesh::MeshDataAccessCollection::TRIANGLES;
            mesh.indices = idx_data;
            mesh.attributes.push_back(vert_attr);
            mesh.attributes.push_back(normal_attr);

            mesh_ = std::make_shared<mesh::MeshDataAccessCollection>();
            mesh_->addMesh("contact_surface", mesh);
        }

        in_data_hash_ = in_call->DataHash();
        frame_id_ = in_call->FrameID();
        ++out_data_hash_;
    }

    out_call->setData(mesh_, out_data_hash_);

    return true;
}


bool megamol::thermodyn::PCLSurface::get_extent_cb(core::Call& c) {
    auto out_call = dynamic_cast<mesh::CallMesh*>(&c);
    if (out_call == nullptr)
        return false;

    auto in_call = in_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_call == nullptr)
        return false;

    auto meta = out_call->getMetaData();

    in_call->SetFrameID(meta.m_frame_ID);
    if (!(*in_call)(1))
        return false;

    meta.m_frame_cnt = in_call->FrameCount();
    meta.m_frame_ID = in_call->FrameID();
    meta.m_bboxs = in_call->AccessBoundingBoxes();

    out_call->setMetaData(meta);

    return true;
}
