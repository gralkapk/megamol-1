#include "PCLSurface.h"

#include "geometry_calls/MultiParticleDataCall.h"
#include "mesh/MeshCalls.h"

#include <pcl/features/normal_3d.h>
#include <pcl/io/pcd_io.h>
#include <pcl/point_types.h>
#include <pcl/search/kdtree.h>
#include <pcl/surface/gp3.h>
#include <pcl/io/obj_io.h>


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

            pcl::GreedyProjectionTriangulation<pcl::PointNormal> gp3;
            pcl::PolygonMesh triangles;

            //gp3.setSearchRadius(0.025);
            gp3.setSearchRadius(6);

            gp3.setMu(2.5);
            gp3.setMaximumNearestNeighbors(10);
            gp3.setMaximumSurfaceAngle(M_PI / 2); // 45 degrees
            gp3.setMinimumAngle(M_PI / 18);       // 10 degrees
            gp3.setMaximumAngle(2 * M_PI / 3);    // 120 degrees
            gp3.setNormalConsistency(false);

            gp3.setInputCloud(cloud_with_normals);
            gp3.setSearchMethod(tree2);
            gp3.reconstruct(triangles);

            std::vector<int> parts = gp3.getPartIDs();
            std::vector<int> states = gp3.getPointStates();

            pcl::io::saveOBJFile("test.obj", triangles);
        }

        in_data_hash_ = in_call->DataHash();
        frame_id_ = in_call->FrameID();
        ++out_data_hash_;
    }

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
