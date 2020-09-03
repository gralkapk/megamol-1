#include "stdafx.h"
#include "ParticleProperties.h"

#include "mmcore/moldyn/MultiParticleDataCall.h"
#include "mmcore/param/EnumParam.h"

#include "adios_plugin/CallADIOSData.h"


megamol::thermodyn::ParticleProperties::ParticleProperties()
    : _data_out_slot("dataOut", "output")
    , _particle_in_slot("particleIn", "input of raw particle data")
    , _density_in_slot("densityIn", "input of density data")
    , _temp_in_slot("tempIn", "input of temperature data")
    , _precision_slot("precision", "set precision of output")
    , _position_order_slot("position_order", "switch between separated and interleaved") {
    _data_out_slot.SetCallback(adios::CallADIOSData::ClassName(), adios::CallADIOSData::FunctionName(0),
        &ParticleProperties::get_data_callback);
    _data_out_slot.SetCallback(adios::CallADIOSData::ClassName(), adios::CallADIOSData::FunctionName(1),
        &ParticleProperties::get_extent_callback);
    MakeSlotAvailable(&_data_out_slot);

    _particle_in_slot.SetCompatibleCall<core::moldyn::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&_particle_in_slot);

    _density_in_slot.SetCompatibleCall<core::moldyn::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&_density_in_slot);

    _temp_in_slot.SetCompatibleCall<core::moldyn::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&_temp_in_slot);

    auto ep = new core::param::EnumParam(static_cast<int>(precision::SINGLE));
    ep->SetTypePair(static_cast<int>(precision::SINGLE), "single");
    ep->SetTypePair(static_cast<int>(precision::DOUBLE), "double");
    _precision_slot << ep;
    MakeSlotAvailable(&_precision_slot);

    ep = new core::param::EnumParam(static_cast<int>(position_order::SEPARATED));
    ep->SetTypePair(static_cast<int>(position_order::SEPARATED), "separated");
    ep->SetTypePair(static_cast<int>(position_order::INTERLEAVED), "interleaved");
    _position_order_slot << ep;
    MakeSlotAvailable(&_position_order_slot);
}


megamol::thermodyn::ParticleProperties::~ParticleProperties() { this->Release(); }


bool megamol::thermodyn::ParticleProperties::create() { return true; }


void megamol::thermodyn::ParticleProperties::release() {}


bool megamol::thermodyn::ParticleProperties::get_data_callback(core::Call& c) {
    auto outCall = dynamic_cast<adios::CallADIOSData*>(&c);
    if (outCall == nullptr) return false;

    core::moldyn::MultiParticleDataCall* particle_in = nullptr;
    core::moldyn::MultiParticleDataCall* density_in = nullptr;
    core::moldyn::MultiParticleDataCall* temp_in = nullptr;

    if (!check_in_calls(particle_in, density_in, temp_in)) return false;

    if (!call_data(particle_in, density_in, temp_in)) return false;

    if (check_recompute(particle_in, density_in, temp_in)) {
        if (!compute_data(particle_in, density_in, temp_in)) return false;
    }

    outCall->setData(_data_map);

    return true;
}


bool megamol::thermodyn::ParticleProperties::get_extent_callback(core::Call& c) {
    auto outCall = dynamic_cast<adios::CallADIOSData*>(&c);
    if (outCall == nullptr) return false;

    core::moldyn::MultiParticleDataCall* particle_in = nullptr;
    core::moldyn::MultiParticleDataCall* density_in = nullptr;
    core::moldyn::MultiParticleDataCall* temp_in = nullptr;

    if (!check_in_calls(particle_in, density_in, temp_in)) return false;

    auto frame_id = outCall->getTime();
    auto force = true;

    if (!call_extents(frame_id, force, particle_in, density_in, temp_in)) return false;


    outCall->setFrameCount(particle_in->FrameCount());
    //outCall->AccessBoundingBoxes() = particle_in->AccessBoundingBoxes();


    return true;
}


bool megamol::thermodyn::ParticleProperties::check_in_calls(core::moldyn::MultiParticleDataCall* particle_in,
    core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in) {
    particle_in = _particle_in_slot.CallAs<core::moldyn::MultiParticleDataCall>();
    if (particle_in == nullptr) return false;

    density_in = _density_in_slot.CallAs<core::moldyn::MultiParticleDataCall>();
    if (density_in == nullptr) {
        megamol::core::utility::log::Log::DefaultLog.WriteWarn("ParticleProperties: Density call not connected");
    }

    temp_in = _temp_in_slot.CallAs<core::moldyn::MultiParticleDataCall>();
    if (temp_in == nullptr) {
        megamol::core::utility::log::Log::DefaultLog.WriteWarn("ParticleProperties: Temperature call not connected");
    }

    return true;
}


bool megamol::thermodyn::ParticleProperties::is_dirty() {
    return _precision_slot.IsDirty() || _position_order_slot.IsDirty();
}


void megamol::thermodyn::ParticleProperties::reset_dirty() {
    _precision_slot.ResetDirty();
    _position_order_slot.ResetDirty();
}


bool megamol::thermodyn::ParticleProperties::call_extents(unsigned int frame_id, bool force,
    core::moldyn::MultiParticleDataCall* particle_in, core::moldyn::MultiParticleDataCall* density_in,
    core::moldyn::MultiParticleDataCall* temp_in) {
    particle_in->SetFrameID(frame_id, force);
    if (!(*particle_in)(1)) return false;

    if (density_in != nullptr) {
        density_in->SetFrameID(frame_id, force);
        if (!(*density_in)(1)) return false;
    }

    if (temp_in != nullptr) {
        temp_in->SetFrameID(frame_id, force);
        if (!(*temp_in)(1)) return false;
    }

    return true;
}


bool megamol::thermodyn::ParticleProperties::call_data(core::moldyn::MultiParticleDataCall* particle_in,
    core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in) {
    if (!(*particle_in)(0)) return false;

    if (density_in != nullptr) {
        if (!(*density_in)(0)) return false;
    }

    if (temp_in != nullptr) {
        if (!(*temp_in)(0)) return false;
    }

    return true;
}


bool megamol::thermodyn::ParticleProperties::compute_data(core::moldyn::MultiParticleDataCall* particle_in,
    core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in) {
    auto density_available = density_in != nullptr;
    auto temp_available = temp_in != nullptr;

    if (particle_in == nullptr) return false;

    // check properties of the inputs

    auto const pl_count = particle_in->GetParticleListCount();

    if (density_available) {
        if (pl_count != density_in->GetParticleListCount()) {
            density_available = false;
            megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                "ParticleProperties: Density is available, but not compatible.");
        }
    }

    if (temp_available) {
        if (pl_count != temp_in->GetParticleListCount()) {
            temp_available = false;
            megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                "ParticleProperties: Temperature is available, but not compatible.");
        }
    }

    auto const prec = static_cast<precision>(_precision_slot.Param<core::param::EnumParam>()->Value());
    auto const order = static_cast<position_order>(_position_order_slot.Param<core::param::EnumParam>()->Value());

    if (prec == precision::SINGLE) {
        std::shared_ptr<adios::FloatContainer> xCont;
        std::shared_ptr<adios::FloatContainer> yCont;
        std::shared_ptr<adios::FloatContainer> zCont;

        std::shared_ptr<adios::FloatContainer> xyzCont;

        if (order == position_order::SEPARATED) {
            xCont = std::make_shared<adios::FloatContainer>();
            yCont = std::make_shared<adios::FloatContainer>();
            zCont = std::make_shared<adios::FloatContainer>();
        } else {
            xyzCont = std::make_shared<adios::FloatContainer>();
        }

        auto pliCont = std::make_shared<adios::UInt32Container>();
        auto radCont = std::make_shared<adios::FloatContainer>();
        auto rCont = std::make_shared<adios::UCharContainer>();
        auto gCont = std::make_shared<adios::UCharContainer>();
        auto bCont = std::make_shared<adios::UCharContainer>();
        auto aCont = std::make_shared<adios::UCharContainer>();
        auto idCont = std::make_shared<adios::UInt64Container>();
        auto dxCont = std::make_shared<adios::FloatContainer>();
        auto dyCont = std::make_shared<adios::FloatContainer>();
        auto dzCont = std::make_shared<adios::FloatContainer>();

        auto densCont = std::make_shared<adios::FloatContainer>();
        auto tempCont = std::make_shared<adios::FloatContainer>();

        std::vector<float> vec_xCont, vec_yCont, vec_zCont, vec_xyzCont;
        if (order == position_order::SEPARATED) {
            vec_xCont = xCont->getVec();
            vec_yCont = yCont->getVec();
            vec_zCont = zCont->getVec();
        } else {
            vec_xyzCont = xyzCont->getVec();
        }

        auto vec_pliCont = pliCont->getVec();
        auto vec_radCont = radCont->getVec();
        auto vec_rCont = rCont->getVec();
        auto vec_gCont = gCont->getVec();
        auto vec_bCont = bCont->getVec();
        auto vec_aCont = aCont->getVec();
        auto vec_idCont = idCont->getVec();
        auto vec_dxCont = dxCont->getVec();
        auto vec_dyCont = dyCont->getVec();
        auto vec_dzCont = dyCont->getVec();

        auto vec_densCont = densCont->getVec();
        auto vec_tempCont = tempCont->getVec();

        for (unsigned int plidx = 0; plidx < pl_count; ++plidx) {
            auto const& plist = particle_in->AccessParticles(plidx);
            auto const p_count = plist.GetCount();

            vec_pliCont.reserve(vec_pliCont.size() + p_count);
            if (order == position_order::SEPARATED) {
                vec_xCont.reserve(vec_xCont.size() + p_count);
                vec_yCont.reserve(vec_yCont.size() + p_count);
                vec_zCont.reserve(vec_zCont.size() + p_count);
            } else {
                vec_xyzCont.reserve(vec_xyzCont.size() + 3 * p_count);
            }
            vec_radCont.reserve(vec_radCont.size() + p_count);
            vec_rCont.reserve(vec_rCont.size() + p_count);
            vec_gCont.reserve(vec_gCont.size() + p_count);
            vec_bCont.reserve(vec_bCont.size() + p_count);
            vec_aCont.reserve(vec_aCont.size() + p_count);
            vec_idCont.reserve(vec_idCont.size() + p_count);
            vec_dxCont.reserve(vec_dxCont.size() + p_count);
            vec_dyCont.reserve(vec_dyCont.size() + p_count);
            vec_dzCont.reserve(vec_dzCont.size() + p_count);

            if (density_available) {
                if (p_count != density_in->AccessParticles(plidx).GetCount()) {
                    density_available = false;
                    megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                        "ParticleProperties: Density is available, but not compatible.");
                } else {
                    vec_densCont.reserve(vec_densCont.size() + p_count);
                }
            }

            if (temp_available) {
                if (p_count != temp_in->AccessParticles(plidx).GetCount()) {
                    temp_available = false;
                    megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                        "ParticleProperties: Temperature is available, but not compatible.");
                } else {
                    vec_tempCont.reserve(vec_tempCont.size() + p_count);
                }
            }

            auto const x_acc = plist.GetParticleStore().GetXAcc();
            auto const y_acc = plist.GetParticleStore().GetYAcc();
            auto const z_acc = plist.GetParticleStore().GetZAcc();
            auto const rad_acc = plist.GetParticleStore().GetRAcc();
            auto const cr_acc = plist.GetParticleStore().GetCRAcc();
            auto const cg_acc = plist.GetParticleStore().GetCGAcc();
            auto const cb_acc = plist.GetParticleStore().GetCBAcc();
            auto const ca_acc = plist.GetParticleStore().GetCAAcc();
            auto const id_acc = plist.GetParticleStore().GetIDAcc();
            auto const dx_acc = plist.GetParticleStore().GetDXAcc();
            auto const dy_acc = plist.GetParticleStore().GetDYAcc();
            auto const dz_acc = plist.GetParticleStore().GetDZAcc();

            if (order == position_order::SEPARATED) {
                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_pliCont.push_back(plidx);
                    vec_xCont.push_back(x_acc->Get_f(pidx));
                    vec_yCont.push_back(y_acc->Get_f(pidx));
                    vec_zCont.push_back(z_acc->Get_f(pidx));
                    vec_radCont.push_back(rad_acc->Get_f(pidx));
                    vec_rCont.push_back(cr_acc->Get_f(pidx));
                    vec_gCont.push_back(cg_acc->Get_f(pidx));
                    vec_bCont.push_back(cb_acc->Get_f(pidx));
                    vec_aCont.push_back(ca_acc->Get_f(pidx));
                    vec_idCont.push_back(id_acc->Get_f(pidx));
                    vec_dxCont.push_back(dx_acc->Get_f(pidx));
                    vec_dyCont.push_back(dy_acc->Get_f(pidx));
                    vec_dzCont.push_back(dz_acc->Get_f(pidx));
                }
            } else {
                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_pliCont.push_back(plidx);
                    vec_xyzCont.push_back(x_acc->Get_f(pidx));
                    vec_xyzCont.push_back(y_acc->Get_f(pidx));
                    vec_xyzCont.push_back(z_acc->Get_f(pidx));
                    vec_radCont.push_back(rad_acc->Get_f(pidx));
                    vec_rCont.push_back(cr_acc->Get_f(pidx));
                    vec_gCont.push_back(cg_acc->Get_f(pidx));
                    vec_bCont.push_back(cb_acc->Get_f(pidx));
                    vec_aCont.push_back(ca_acc->Get_f(pidx));
                    vec_idCont.push_back(id_acc->Get_f(pidx));
                    vec_dxCont.push_back(dx_acc->Get_f(pidx));
                    vec_dyCont.push_back(dy_acc->Get_f(pidx));
                    vec_dzCont.push_back(dz_acc->Get_f(pidx));
                }
            }

            if (density_available) {
                auto const dens_acc = density_in->AccessParticles(plidx).GetParticleStore().GetCRAcc();

                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_densCont.push_back(dens_acc->Get_f(pidx));
                }
            }

            if (temp_available) {
                auto const temp_acc = temp_in->AccessParticles(plidx).GetParticleStore().GetCRAcc();

                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_tempCont.push_back(temp_acc->Get_f(pidx));
                }
            }

            _data_map->clear();

            (*_data_map)["plidx"] = std::move(pliCont);
            if (order == position_order::SEPARATED) {
                (*_data_map)["x"] = std::move(xCont);
                (*_data_map)["y"] = std::move(yCont);
                (*_data_map)["z"] = std::move(zCont);
            } else {
                (*_data_map)["xyz"] = std::move(xyzCont);
            }
            (*_data_map)["rad"] = std::move(radCont);
            (*_data_map)["cr"] = std::move(rCont);
            (*_data_map)["cg"] = std::move(gCont);
            (*_data_map)["cb"] = std::move(bCont);
            (*_data_map)["ca"] = std::move(aCont);
            (*_data_map)["id"] = std::move(idCont);
            if (plist.GetDirDataType() != core::moldyn::SimpleSphericalParticles::DIRDATA_NONE) {
                (*_data_map)["dr"] = std::move(dxCont);
                (*_data_map)["dg"] = std::move(dyCont);
                (*_data_map)["db"] = std::move(dzCont);
            }
            if (density_available) {
                (*_data_map)["dens"] = std::move(densCont);
            }
            if (temp_available) {
                (*_data_map)["temp"] = std::move(tempCont);
            }
        }
    } else {
        std::shared_ptr<adios::DoubleContainer> xCont;
        std::shared_ptr<adios::DoubleContainer> yCont;
        std::shared_ptr<adios::DoubleContainer> zCont;

        std::shared_ptr<adios::DoubleContainer> xyzCont;

        if (order == position_order::SEPARATED) {
            xCont = std::make_shared<adios::DoubleContainer>();
            yCont = std::make_shared<adios::DoubleContainer>();
            zCont = std::make_shared<adios::DoubleContainer>();
        } else {
            xyzCont = std::make_shared<adios::DoubleContainer>();
        }

        auto pliCont = std::make_shared<adios::UInt32Container>();
        auto radCont = std::make_shared<adios::DoubleContainer>();
        auto rCont = std::make_shared<adios::UCharContainer>();
        auto gCont = std::make_shared<adios::UCharContainer>();
        auto bCont = std::make_shared<adios::UCharContainer>();
        auto aCont = std::make_shared<adios::UCharContainer>();
        auto idCont = std::make_shared<adios::UInt64Container>();
        auto dxCont = std::make_shared<adios::DoubleContainer>();
        auto dyCont = std::make_shared<adios::DoubleContainer>();
        auto dzCont = std::make_shared<adios::DoubleContainer>();

        auto densCont = std::make_shared<adios::DoubleContainer>();
        auto tempCont = std::make_shared<adios::DoubleContainer>();

        std::vector<double> vec_xCont, vec_yCont, vec_zCont, vec_xyzCont;
        if (order == position_order::SEPARATED) {
            vec_xCont = xCont->getVec();
            vec_yCont = yCont->getVec();
            vec_zCont = zCont->getVec();
        } else {
            vec_xyzCont = xyzCont->getVec();
        }

        auto vec_pliCont = pliCont->getVec();
        auto vec_radCont = radCont->getVec();
        auto vec_rCont = rCont->getVec();
        auto vec_gCont = gCont->getVec();
        auto vec_bCont = bCont->getVec();
        auto vec_aCont = aCont->getVec();
        auto vec_idCont = idCont->getVec();
        auto vec_dxCont = dxCont->getVec();
        auto vec_dyCont = dyCont->getVec();
        auto vec_dzCont = dyCont->getVec();

        auto vec_densCont = densCont->getVec();
        auto vec_tempCont = tempCont->getVec();

        for (unsigned int plidx = 0; plidx < pl_count; ++plidx) {
            auto const& plist = particle_in->AccessParticles(plidx);
            auto const p_count = plist.GetCount();

            vec_pliCont.reserve(vec_pliCont.size() + p_count);
            if (order == position_order::SEPARATED) {
                vec_xCont.reserve(vec_xCont.size() + p_count);
                vec_yCont.reserve(vec_yCont.size() + p_count);
                vec_zCont.reserve(vec_zCont.size() + p_count);
            } else {
                vec_xyzCont.reserve(vec_xyzCont.size() + 3 * p_count);
            }
            vec_radCont.reserve(vec_radCont.size() + p_count);
            vec_rCont.reserve(vec_rCont.size() + p_count);
            vec_gCont.reserve(vec_gCont.size() + p_count);
            vec_bCont.reserve(vec_bCont.size() + p_count);
            vec_aCont.reserve(vec_aCont.size() + p_count);
            vec_idCont.reserve(vec_idCont.size() + p_count);
            vec_dxCont.reserve(vec_dxCont.size() + p_count);
            vec_dyCont.reserve(vec_dyCont.size() + p_count);
            vec_dzCont.reserve(vec_dzCont.size() + p_count);

            if (density_available) {
                if (p_count != density_in->AccessParticles(plidx).GetCount()) {
                    density_available = false;
                    megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                        "ParticleProperties: Density is available, but not compatible.");
                } else {
                    vec_densCont.reserve(vec_densCont.size() + p_count);
                }
            }

            if (temp_available) {
                if (p_count != temp_in->AccessParticles(plidx).GetCount()) {
                    temp_available = false;
                    megamol::core::utility::log::Log::DefaultLog.WriteWarn(
                        "ParticleProperties: Temperature is available, but not compatible.");
                } else {
                    vec_tempCont.reserve(vec_tempCont.size() + p_count);
                }
            }

            auto const x_acc = plist.GetParticleStore().GetXAcc();
            auto const y_acc = plist.GetParticleStore().GetYAcc();
            auto const z_acc = plist.GetParticleStore().GetZAcc();
            auto const rad_acc = plist.GetParticleStore().GetRAcc();
            auto const cr_acc = plist.GetParticleStore().GetCRAcc();
            auto const cg_acc = plist.GetParticleStore().GetCGAcc();
            auto const cb_acc = plist.GetParticleStore().GetCBAcc();
            auto const ca_acc = plist.GetParticleStore().GetCAAcc();
            auto const id_acc = plist.GetParticleStore().GetIDAcc();
            auto const dx_acc = plist.GetParticleStore().GetDXAcc();
            auto const dy_acc = plist.GetParticleStore().GetDYAcc();
            auto const dz_acc = plist.GetParticleStore().GetDZAcc();

            if (order == position_order::SEPARATED) {
                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_pliCont.push_back(plidx);
                    vec_xCont.push_back(x_acc->Get_d(pidx));
                    vec_yCont.push_back(y_acc->Get_d(pidx));
                    vec_zCont.push_back(z_acc->Get_d(pidx));
                    vec_radCont.push_back(rad_acc->Get_d(pidx));
                    vec_rCont.push_back(cr_acc->Get_u8(pidx));
                    vec_gCont.push_back(cg_acc->Get_u8(pidx));
                    vec_bCont.push_back(cb_acc->Get_u8(pidx));
                    vec_aCont.push_back(ca_acc->Get_u8(pidx));
                    vec_idCont.push_back(id_acc->Get_u64(pidx));
                    vec_dxCont.push_back(dx_acc->Get_d(pidx));
                    vec_dyCont.push_back(dy_acc->Get_d(pidx));
                    vec_dzCont.push_back(dz_acc->Get_d(pidx));
                }
            } else {
                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_pliCont.push_back(plidx);
                    vec_xyzCont.push_back(x_acc->Get_d(pidx));
                    vec_xyzCont.push_back(y_acc->Get_d(pidx));
                    vec_xyzCont.push_back(z_acc->Get_d(pidx));
                    vec_radCont.push_back(rad_acc->Get_d(pidx));
                    vec_rCont.push_back(cr_acc->Get_u8(pidx));
                    vec_gCont.push_back(cg_acc->Get_u8(pidx));
                    vec_bCont.push_back(cb_acc->Get_u8(pidx));
                    vec_aCont.push_back(ca_acc->Get_u8(pidx));
                    vec_idCont.push_back(id_acc->Get_u64(pidx));
                    vec_dxCont.push_back(dx_acc->Get_d(pidx));
                    vec_dyCont.push_back(dy_acc->Get_d(pidx));
                    vec_dzCont.push_back(dz_acc->Get_d(pidx));
                }
            }

            if (density_available) {
                auto const dens_acc = density_in->AccessParticles(plidx).GetParticleStore().GetCRAcc();

                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_densCont.push_back(dens_acc->Get_d(pidx));
                }
            }

            if (temp_available) {
                auto const temp_acc = temp_in->AccessParticles(plidx).GetParticleStore().GetCRAcc();

                for (size_t pidx = 0; pidx < p_count; ++pidx) {
                    vec_tempCont.push_back(temp_acc->Get_d(pidx));
                }
            }

            _data_map->clear();

            (*_data_map)["plidx"] = std::move(pliCont);
            if (order == position_order::SEPARATED) {
                (*_data_map)["x"] = std::move(xCont);
                (*_data_map)["y"] = std::move(yCont);
                (*_data_map)["z"] = std::move(zCont);
            } else {
                (*_data_map)["xyz"] = std::move(xyzCont);
            }
            (*_data_map)["rad"] = std::move(radCont);
            (*_data_map)["cr"] = std::move(rCont);
            (*_data_map)["cg"] = std::move(gCont);
            (*_data_map)["cb"] = std::move(bCont);
            (*_data_map)["ca"] = std::move(aCont);
            (*_data_map)["id"] = std::move(idCont);
            if (plist.GetDirDataType() != core::moldyn::SimpleSphericalParticles::DIRDATA_NONE) {
                (*_data_map)["dr"] = std::move(dxCont);
                (*_data_map)["dg"] = std::move(dyCont);
                (*_data_map)["db"] = std::move(dzCont);
            }
            if (density_available) {
                (*_data_map)["dens"] = std::move(densCont);
            }
            if (temp_available) {
                (*_data_map)["temp"] = std::move(tempCont);
            }
        }
    }
}


bool megamol::thermodyn::ParticleProperties::check_recompute(core::moldyn::MultiParticleDataCall* particle_in,
    core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in) {
    bool recompute = false;
    if (is_dirty()) recompute = true;

    if (particle_in->DataHash() != _particle_in_data_hash || particle_in->FrameID() != _in_frame_id) {
        _particle_in_data_hash = particle_in->DataHash();
        recompute = true;
    }

    if (density_in != nullptr) {
        if (density_in->DataHash() != _density_in_data_hash || density_in->FrameID() != _in_frame_id) {
            _density_in_data_hash = density_in->DataHash();
            recompute = true;
        }
    }

    if (temp_in != nullptr) {
        if (temp_in->DataHash() != _density_in_data_hash || temp_in->FrameID() != _in_frame_id) {
            _temp_in_data_hash = temp_in->DataHash();
            recompute = true;
        }
    }

    return recompute;
}
