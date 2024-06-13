#include "PKDUtils.h"

#include <glm/glm.hpp>
#include <tbb/parallel_for.h>

#include "pkd_utils.h"

namespace megamol::optix_hpg {
// BEGIN PKD



void recBuild(size_t /* root node */ P, device::PKDParticle* particle, size_t N, device::box3f bounds, device::FPKDParticle* pack = nullptr) {
    if (P >= N)
        return;

    int dim = arg_max(bounds.upper - bounds.lower);

    const size_t L = lChild(P);
    const size_t R = rChild(P);
    const bool lValid = (L < N);
    const bool rValid = (R < N);
    makeHeap(std::greater<float>(), L, particle, N, dim, pack);
    makeHeap(std::less<float>(), R, particle, N, dim, pack);

    if (rValid) {
        while (particle[L].pos[dim] > particle[R].pos[dim]) {
            std::swap(particle[L], particle[R]);
            if (pack) {
                std::swap(pack[L], pack[R]);
            }
            trickle(std::greater<float>(), L, particle, N, dim, pack);
            trickle(std::less<float>(), R, particle, N, dim, pack);
        }
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            if (pack) {
                std::swap(pack[L], pack[P]);
                pack[L].dim = dim;
            }
            particle[L].dim = dim;
        } else if (particle[R].pos[dim] < particle[P].pos[dim]) {
            std::swap(particle[R], particle[P]);
            if (pack) {
                std::swap(pack[R], pack[P]);
                pack[R].dim = dim;
            }
            particle[R].dim = dim;
        } else
            /* nothing, root fits */;
    } else if (lValid) {
        if (particle[L].pos[dim] > particle[P].pos[dim]) {
            std::swap(particle[L], particle[P]);
            if (pack) {
                std::swap(pack[L], pack[P]);
                pack[L].dim = dim;
            }
            particle[L].dim = dim;
        }
    }

    device::box3f lBounds = bounds;
    device::box3f rBounds = bounds;
    lBounds.upper[dim] = rBounds.lower[dim] = particle[P].pos[dim];
    particle[P].dim = dim;

    tbb::parallel_for(0, 2, [&](int childID) {
        if (childID) {
            recBuild(L, particle, N, lBounds, pack);
        } else {
            recBuild(R, particle, N, rBounds, pack);
        }
    });
}

void makePKD(std::vector<device::PKDParticle>& particles, device::box3f bounds) {
    recBuild(/*node:*/ 0, particles.data(), particles.size(), bounds);
}

void makePKD(std::vector<device::PKDParticle>& particles, size_t begin, size_t end, device::box3f bounds, device::FPKDParticle* pack) {
    recBuild(/*node:*/ 0, particles.data() + begin, end - begin, bounds, pack);
}



void recBuild(size_t /* root node */ P, device::SPKDParticle* particle, size_t N, device::box3f const& bounds, device::SPKDlet const& treelet) {
    if (P >= N)
        return;

    int dim = arg_max(bounds.upper - bounds.lower);

    const size_t L = lChild(P);
    const size_t R = rChild(P);
    const bool lValid = (L < N);
    const bool rValid = (R < N);
    makeHeap(std::greater<float>(), L, particle, N, dim, treelet);
    makeHeap(std::less<float>(), R, particle, N, dim, treelet);

    auto const P_pos = decode_spart(particle[P], treelet);

    if (rValid) {
        while (decode_spart(particle[L], treelet)[dim] > decode_spart(particle[R], treelet)[dim]) {
            std::swap(particle[L], particle[R]);
            trickle(std::greater<float>(), L, particle, N, dim, treelet);
            trickle(std::less<float>(), R, particle, N, dim, treelet);
        }
        if (decode_spart(particle[L], treelet)[dim] > P_pos[dim]) {
            std::swap(particle[L], particle[P]);
            particle[L].dim = dim;
        } else if (decode_spart(particle[R], treelet)[dim] < P_pos[dim]) {
            std::swap(particle[R], particle[P]);
            particle[R].dim = dim;
        } else
            /* nothing, root fits */;
    } else if (lValid) {
        if (decode_spart(particle[L], treelet)[dim] > P_pos[dim]) {
            std::swap(particle[L], particle[P]);
            particle[L].dim = dim;
        }
    }

    device::box3f lBounds = bounds;
    device::box3f rBounds = bounds;
    lBounds.upper[dim] = rBounds.lower[dim] = P_pos[dim];
    particle[P].dim = dim;

    tbb::parallel_for(0, 2, [&](int childID) {
        if (childID) {
            recBuild(L, particle, N, lBounds, treelet);
        } else {
            recBuild(R, particle, N, rBounds, treelet);
        }
    });
}

void makePKD(std::vector<device::SPKDParticle>& particles, device::SPKDlet const& treelet, size_t begin) {
    recBuild(/*node:*/ 0, particles.data() + treelet.begin - begin, treelet.end - treelet.begin, treelet.bounds, treelet);
}

// END PKD


// BEGIN TREELETS

//size_t sort_partition(
//    std::vector<device::PKDParticle>& particles, size_t begin, size_t end, device::box3f bounds, int& splitDim) {
//    // -------------------------------------------------------
//    // determine split pos
//    // -------------------------------------------------------
//    splitDim = arg_max(bounds.upper - bounds.lower);
//    //float splitPos = bounds.center()[splitDim];
//    float splitPos = (0.5f * (bounds.upper + bounds.lower))[splitDim];
//
//    // -------------------------------------------------------
//    // now partition ...
//    // -------------------------------------------------------
//    size_t mid = begin;
//    size_t l = begin, r = (end - 1);
//    // quicksort partition:
//    while (l <= r) {
//        while (l < r && particles[l].pos[splitDim] < splitPos)
//            ++l;
//        while (l < r && particles[r].pos[splitDim] >= splitPos)
//            --r;
//        if (l == r) {
//            mid = l;
//            break;
//        }
//
//        std::swap(particles[l], particles[r]);
//    }
//
//    // catch-all for extreme cases where all particles are on the same
//    // spot, and can't be split:
//    if (mid == begin || mid == end)
//        mid = (begin + end) / 2;
//
//    return mid;
//}

device::box3f extendBounds(std::vector<device::PKDParticle> const& particles, size_t begin, size_t end, float radius) {
    device::box3f bounds;
    for (int64_t p_idx = begin; p_idx < end; ++p_idx) {
        auto const new_lower = particles[p_idx].pos - radius;
        auto const new_upper = particles[p_idx].pos + radius;
        bounds.extend(new_lower);
        bounds.extend(new_upper);
    }

    return bounds;
}

std::vector<device::PKDlet> prePartition_inPlace(std::vector<device::PKDParticle>& particles, size_t maxSize,
    float radius, std::function<bool(device::box3f const&)> add_cond) {
    std::mutex resultMutex;
    std::vector<device::PKDlet> result;

    if (add_cond == nullptr) {
        partitionRecursively<device::PKDParticle, device::box3f>(
            particles, 0ULL, particles.size(), [&](size_t begin, size_t end, bool force, device::box3f const& bounds) {
                /*bool makeLeaf() :*/
                const size_t size = end - begin;
                if (size > maxSize && !force)
                    return false;

                device::PKDlet treelet;
                treelet.begin = begin;
                treelet.end = end;
                //treelet.bounds = extendBounds(particles, begin, end, radius);
                treelet.bounds = bounds;
                treelet.bounds.lower -= radius;
                treelet.bounds.upper += radius;

                std::lock_guard<std::mutex> lock(resultMutex);
                result.push_back(treelet);
                return true;
            });
    } else {
        partitionRecursively<device::PKDParticle, device::box3f>(
            particles, 0ULL, particles.size(), [&](size_t begin, size_t end, bool force, device::box3f const& bounds) {
                /*bool makeLeaf() :*/
                const size_t size = end - begin;
                if ((size > maxSize && !force) || add_cond(bounds))
                    return false;

                device::PKDlet treelet;
                treelet.begin = begin;
                treelet.end = end;
                //treelet.bounds = extendBounds(particles, begin, end, radius);
                treelet.bounds = bounds;
                treelet.bounds.lower -= radius;
                treelet.bounds.upper += radius;

                std::lock_guard<std::mutex> lock(resultMutex);
                result.push_back(treelet);
                return true;
            });
    }

    return std::move(result);
}

// END TREELETS

// BEGIN COMPRESSION

//glm::uvec3 encode_coord(glm::vec3 const& pos, glm::vec3 const& center, glm::vec3 const& span) {
//    constexpr unsigned int digits = 1023u;
//    auto const dir = pos - center;
//    auto const diff = span / static_cast<float>(digits);
//    auto const coord = glm::uvec3(dir.x / diff.x, dir.y / diff.y, dir.z / diff.z);
//
//    return coord;
//}
//
//glm::vec3 decode_coord(glm::uvec3 const& coord, glm::vec3 const& center, glm::vec3 const& span) {
//    constexpr unsigned int digits = 1023u;
//    auto const diff = span / static_cast<float>(digits);
//    auto pos = glm::vec3(static_cast<float>(coord.x) * diff.x, static_cast<float>(coord.y) * diff.y,
//        static_cast<float>(coord.z) * diff.z);
//    pos = pos + center;
//
//    return pos;
//}

std::tuple<std::vector<device::PKDlet>, std::vector<std::pair<unsigned int, device::QPKDParticle>>> makeSpartition(
    std::vector<device::QPKDParticle> const& data, size_t begin, size_t end, float radius) {
    std::vector<std::pair<unsigned int, device::QPKDParticle>> particles(end - begin);
    for (size_t i = begin; i < end; ++i) {
        particles[i - begin] = std::make_pair(i, data[i]);
    }
    std::vector<device::PKDlet> treelets;
    treelets.reserve(particles.size() / 256);
    // divide in x axis
    std::sort(particles.begin(), particles.end(),
        [](auto const& lhs, auto const& rhs) { return lhs.second.x < rhs.second.x; });
    auto sec_it = particles.begin();
    while (sec_it != particles.end()) {
        byte_cast fbc;
        fbc.ui = sec_it->second.x;
        auto sec_begin = sec_it;
        auto sec_end = std::find_if_not(sec_begin, particles.end(), [&fbc](auto const& el) {
            byte_cast bc;
            bc.ui = el.second.x;
            return fbc.parts.b == bc.parts.b;
        });
        // we got the range of equal prefix
        device::PKDlet treelet;
        treelet.begin = std::distance(particles.begin(), sec_begin);
        treelet.end = std::distance(particles.begin(), sec_end);
        treelets.push_back(treelet);
        sec_it = sec_end;
    }
    // divide in y axis
    std::vector<device::PKDlet> tmp_treelets;
    tmp_treelets.reserve(particles.size() / 256);
    for (auto const& el : treelets) {
        std::sort(particles.begin() + el.begin, particles.begin() + el.end,
            [](auto const& lhs, auto const& rhs) { return lhs.second.y < rhs.second.y; });
        auto sec_it = particles.begin() + el.begin;
        while (sec_it != (particles.begin() + el.end)) {
            byte_cast fbc;
            fbc.ui = sec_it->second.y;
            auto sec_begin = sec_it;
            auto sec_end = std::find_if_not(sec_begin, particles.begin() + el.end, [&fbc](auto const& e) {
                byte_cast bc;
                bc.ui = e.second.y;
                return fbc.parts.b == bc.parts.b;
            });
            device::PKDlet treelet;
            treelet.begin = std::distance(particles.begin(), sec_begin);
            treelet.end = std::distance(particles.begin(), sec_end);
            tmp_treelets.push_back(treelet);
            sec_it = sec_end;
        }
    }
    treelets = tmp_treelets;
    tmp_treelets.clear();
    tmp_treelets.reserve(particles.size() / 256);
    // divide in z axis
    for (auto const& el : treelets) {
        std::sort(particles.begin() + el.begin, particles.begin() + el.end,
            [](auto const& lhs, auto const& rhs) { return lhs.second.z < rhs.second.z; });
        auto sec_it = particles.begin() + el.begin;
        while (sec_it != (particles.begin() + el.end)) {
            byte_cast fbc;
            fbc.ui = sec_it->second.z;
            auto sec_begin = sec_it;
            auto sec_end = std::find_if_not(sec_begin, particles.begin() + el.end, [&fbc](auto const& e) {
                byte_cast bc;
                bc.ui = e.second.z;
                return fbc.parts.b == bc.parts.b;
            });
            device::PKDlet treelet;
            treelet.begin = std::distance(particles.begin(), sec_begin);
            treelet.end = std::distance(particles.begin(), sec_end);
            tmp_treelets.push_back(treelet);
            sec_it = sec_end;
        }
    }
    treelets = tmp_treelets;
    for (auto& el : treelets) {
        el.begin += begin;
        el.end += begin;
    }
    return std::make_tuple(treelets, particles);
}

//std::tuple<std::vector<device::SPKDlet>, std::vector<device::SPKDParticle>> slice_qparticles(
//    std::vector<device::PKDlet> const& treelets,
//    std::vector<std::pair<unsigned int, device::QPKDParticle>> const& particles,
//    std::vector<device::PKDParticle> const& org_data, size_t begin, size_t end, float radius) {
//    std::vector<device::SPKDlet> streelets;
//    streelets.reserve(treelets.size());
//    std::vector<device::SPKDParticle> sparticles(particles.size());
//    for (auto const& treelet : treelets) {
//        std::vector<device::PKDParticle> tmp_data(treelet.end - treelet.begin);
//        for (size_t i = treelet.begin; i < treelet.end; ++i) {
//            //tmp_data[i - treelet.begin].pos = decode_coord(particles[i], glm::vec3(), glm::vec3());
//            tmp_data[i - treelet.begin].pos = org_data[particles[i - begin].first].pos;
//            byte_cast bcx;
//            bcx.ui = particles[i - begin].second.x;
//            sparticles[i - begin].x = bcx.parts.a;
//            byte_cast bcy;
//            bcy.ui = particles[i - begin].second.y;
//            sparticles[i - begin].y = bcy.parts.a;
//            byte_cast bcz;
//            bcz.ui = particles[i - begin].second.z;
//            sparticles[i - begin].z = bcz.parts.a;
//        }
//        auto const bounds = extendBounds(tmp_data, 0, tmp_data.size(), radius);
//        device::SPKDlet st;
//        st.begin = treelet.begin;
//        st.end = treelet.end;
//        st.bounds = bounds;
//        byte_cast bcx;
//        bcx.ui = particles[treelet.begin - begin].second.x;
//        st.sx = bcx.parts.b;
//        byte_cast bcy;
//        bcy.ui = particles[treelet.begin - begin].second.y;
//        st.sy = bcy.parts.b;
//        byte_cast bcz;
//        bcz.ui = particles[treelet.begin - begin].second.z;
//        st.sz = bcz.parts.b;
//        streelets.push_back(st);
//    }
//    return std::make_tuple(streelets, sparticles);
//}
//
//std::vector<glm::vec3> compute_diffs(std::vector<device::SPKDlet> const& treelets,
//    std::vector<device::SPKDParticle> const& sparticles,
//    std::vector<std::pair<unsigned int, device::QPKDParticle>> const& qparticles,
//    std::vector<device::PKDParticle> const& org_data, size_t begin, size_t end, glm::vec3 const& lower) {
//    std::vector<glm::vec3> diffs(sparticles.size());
//    for (auto const& treelet : treelets) {
//        for (size_t i = treelet.begin; i < treelet.end; ++i) {
//            device::QPKDParticle qp;
//            byte_cast bc;
//            bc.ui = 0;
//            bc.parts.a = sparticles[i - begin].x;
//            bc.parts.b = treelet.sx;
//            qp.x = bc.ui;
//            bc.parts.a = sparticles[i - begin].y;
//            bc.parts.b = treelet.sy;
//            qp.y = bc.ui;
//            bc.parts.a = sparticles[i - begin].z;
//            bc.parts.b = treelet.sz;
//            qp.z = bc.ui;
//            glm::dvec3 pos = decode_coord(qp /*, glm::vec3(), glm::vec3()*/) + lower;
//            glm::dvec3 qpos = decode_coord(qparticles[i - begin].second /*, glm::vec3(), glm::vec3()*/) + lower;
//            glm::dvec3 org_pos = org_data[qparticles[i - begin].first].pos;
//            diffs[i - begin] = pos - org_pos;
//        }
//    }
//    return diffs;
//}


std::vector<std::pair<size_t, size_t>> gridify(
    std::vector<device::PKDParticle>& data, glm::vec3 const& lower, glm::vec3 const& upper) {
    constexpr float const split_size = (1 << (16 - dec_val));//    -1.0f;
    auto const span = upper - lower;
    auto const num_cells = glm::ceil(span / split_size);
    auto const diff = span / num_cells;
    std::vector<int> cell_idxs(data.size());
    std::vector<size_t> num_elements(num_cells.x * num_cells.y * num_cells.z, 0);
    for (size_t i = 0; i < data.size(); ++i) {
        glm::ivec3 cell_idx = data[i].pos / diff;
        cell_idx.x = cell_idx.x >= num_cells.x ? num_cells.x - 1 : cell_idx.x;
        cell_idx.y = cell_idx.y >= num_cells.y ? num_cells.y - 1 : cell_idx.y;
        cell_idx.z = cell_idx.z >= num_cells.z ? num_cells.z - 1 : cell_idx.z;
        cell_idxs[i] = cell_idx.x + num_cells.x * (cell_idx.y + num_cells.y * cell_idx.z);
        ++num_elements[cell_idxs[i]];
    }
    std::vector<std::pair<size_t, size_t>> grid_cells(num_elements.size(), std::make_pair(0, 0));
    grid_cells[0].second = num_elements[0];
    for (size_t i = 1; i < num_elements.size(); ++i) {
        num_elements[i] += num_elements[i - 1];
        grid_cells[i].first = num_elements[i - 1];
        grid_cells[i].second = num_elements[i];
    }
    std::vector<device::PKDParticle> tmp(data.size());
    for (size_t i = 0; i < data.size(); ++i) {
        auto const idx = cell_idxs[i];
        tmp[--num_elements[idx]] = data[i];
    }
    data = tmp;
    /*auto max_x = (*std::max_element(
        cell_idx.begin(), cell_idx.end(), [](auto const& lhs, auto const& rhs) { return lhs.x < rhs.x; })).x;
    auto max_y = (*std::max_element(
        cell_idx.begin(), cell_idx.end(), [](auto const& lhs, auto const& rhs) { return lhs.y < rhs.y; })).y;
    auto max_z = (*std::max_element(
        cell_idx.begin(), cell_idx.end(), [](auto const& lhs, auto const& rhs) { return lhs.z < rhs.z; })).z;*/
    return grid_cells;
}

void convert(size_t P, device::PKDParticle* in_particle, device::QPKDParticle* out_particle, size_t N,
    device::box3f bounds, float radius, device::PKDParticle* out_decode, glm::uvec3* out_coord) {
    if (P >= N)
        return;

    constexpr float target_prec = 1.0e-5f;

    auto const center = PKD_BOUNDS_CENTER;
    auto const span = bounds.span();
    /*auto dig = span / target_prec;
    dig = glm::log2(dig);*/

    auto coord = encode_coord(in_particle[P].pos, center, span);
    auto const pos = decode_coord(coord /*, center, span*/);
    /*if (out_decode) {
        auto d = glm::dvec3(in_particle[P].pos) - glm::dvec3(pos);
        out_decode[P].pos = d;
    }*/
    if (out_coord) {
        /*constexpr unsigned int digits = 32767u;
        auto dir = pos - center;
        auto const diff = span / static_cast<float>(digits);
        auto const x = static_cast<unsigned int>(dir.x / diff.x);
        auto const y = static_cast<unsigned int>(dir.y / diff.y);
        auto const z = static_cast<unsigned int>(dir.z / diff.z);
        out_coord[P] = glm::uvec3(x, y, z);*/
        //auto dir = in_particle[P].pos - in_particle[parent(P)].pos;
        ////auto const diff = glm::log2(dir / target_prec);
        //auto const diff = glm::uvec3(glm::abs(dir / target_prec));
        //out_coord[P] = glm::uvec3(diff.x, diff.y, diff.z);

        decvec3 dec_pos = in_particle[P].pos;
        decvec3 dec_parent = in_particle[parent(P)].pos;
        decvec3 dec_span = span;
        auto const diff = (dec_pos - dec_parent);
        out_coord[P] = glm::uvec3(diff.x, diff.y, diff.z);
        if (out_decode) {
            glm::vec3 tmp = (diff + dec_parent);
            auto const d = glm::dvec3(in_particle[P].pos) - glm::dvec3(tmp);
            out_decode[P].pos = d;
        }
    }

    int const dim = in_particle[P].dim;
    encode_dim(dim, coord);
    out_particle[P] = coord;
    //out_particle[P].dim = dim;

    auto lBounds = bounds;
    auto rBounds = bounds;

    lBounds.upper[dim] = pos[dim] + radius;
    rBounds.lower[dim] = pos[dim] - radius;

    auto const L = lChild(P);
    auto const R = rChild(P);
    //const bool lValid = (L < N);
    //const bool rValid = (R < N);

    // TODO parallel
    convert(L, in_particle, out_particle, N, bounds, radius, out_decode, out_coord);
    convert(R, in_particle, out_particle, N, bounds, radius, out_decode, out_coord);
}

// END COMPRESSION
} // namespace megamol::optix_hpg
