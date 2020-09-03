#pragma once

#include "mmcore/Call.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

#include "mmcore/moldyn/MultiParticleDataCall.h"

#include "adios_plugin/CallADIOSData.h"

namespace megamol::thermodyn {

class ParticleProperties : public core::Module {
public:
    enum class precision : unsigned char { SINGLE = 0, DOUBLE };

    enum class position_order : unsigned char { SEPARATED = 0, INTERLEAVED };

    /** Return module class name */
    static const char* ClassName(void) { return "ParticleProperties"; }

    /** Return module class description */
    static const char* Description(void) { return "Collects particle properties for output in ADIOS"; }

    /** Module is always available */
    static bool IsAvailable(void) { return true; }

    ParticleProperties();

    virtual ~ParticleProperties();

protected:
    bool create() override;

    void release() override;

private:
    bool get_data_callback(core::Call& c);

    bool get_extent_callback(core::Call& c);

    bool check_in_calls(core::moldyn::MultiParticleDataCall* particle_in,
        core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in);

    bool is_dirty();

    void reset_dirty();

    bool call_extents(unsigned int frame_id, bool force, core::moldyn::MultiParticleDataCall* particle_in,
        core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in);

    bool call_data(core::moldyn::MultiParticleDataCall* particle_in, core::moldyn::MultiParticleDataCall* density_in,
        core::moldyn::MultiParticleDataCall* temp_in);

    bool compute_data(core::moldyn::MultiParticleDataCall* particle_in, core::moldyn::MultiParticleDataCall* density_in,
        core::moldyn::MultiParticleDataCall* temp_in);

    bool check_recompute(core::moldyn::MultiParticleDataCall* particle_in,
        core::moldyn::MultiParticleDataCall* density_in, core::moldyn::MultiParticleDataCall* temp_in);

    core::CalleeSlot _data_out_slot;

    core::CallerSlot _particle_in_slot;

    core::CallerSlot _density_in_slot;

    core::CallerSlot _temp_in_slot;

    core::param::ParamSlot _precision_slot;

    core::param::ParamSlot _position_order_slot;

    size_t _out_data_hash;

    size_t _particle_in_data_hash;

    size_t _density_in_data_hash;

    size_t _temp_in_data_hash;

    unsigned int _out_frame_id;

    unsigned int _in_frame_id;

    std::shared_ptr<adios::adiosDataMap> _data_map;


}; // end class ParticleProperties

} // end namespace megamol::thermodyn
