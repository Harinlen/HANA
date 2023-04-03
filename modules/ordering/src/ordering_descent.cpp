#include <cassert>
#include <limits>
#include <random>
#include <cmath>

#include "hmr_global.hpp"
#include "hmr_ga.hpp"
#include "hmr_ui.hpp"

#include "ordering_descent.hpp"

constexpr double LIMIT = 10000000;
const double LimitLog = log(LIMIT);

typedef struct ORDERING_TIG
{
    int32_t idx;
    int32_t size;
} ORDERING_TIG;

typedef struct ORDERING_GENOME
{
    ORDERING_TIG *tigs;
    size_t tig_size;
    ORDERING_COUNTS *edges;
} ORDERING_GENOME;

void *ordering_genome_create(std::mt19937_64 &rng, void *user)
{
    HMR_UNUSED(rng)
    ORDERING_INFO *info = reinterpret_cast<ORDERING_INFO *>(user);
    ORDERING_GENOME *genome = static_cast<ORDERING_GENOME *>(malloc(sizeof(ORDERING_GENOME)));
    genome->edges = &info->edges;
    genome->tig_size = info->contig_size;
    size_t tig_bytes = sizeof(ORDERING_TIG) * genome->tig_size;
    genome->tigs = static_cast<ORDERING_TIG *>(malloc(tig_bytes));
    memcpy(genome->tigs, info->init_genome, tig_bytes);
    return genome;
}

void ordering_genome_free(void *user)
{
    ORDERING_GENOME *genome = reinterpret_cast<ORDERING_GENOME *>(user);
    free(genome->tigs);
    free(genome);
}

void *ordering_genome_clone(void *src_genome)
{
    ORDERING_GENOME *src = reinterpret_cast<ORDERING_GENOME *>(src_genome);
    ORDERING_GENOME *dst = static_cast<ORDERING_GENOME *>(malloc(sizeof(ORDERING_GENOME)));
    dst->edges = src->edges;
    dst->tig_size = src->tig_size;
    size_t tig_bytes = sizeof(ORDERING_TIG) * dst->tig_size;
    dst->tigs = static_cast<ORDERING_TIG *>(malloc(tig_bytes));
    memcpy(dst->tigs, src->tigs, tig_bytes);
    return dst;
}

void ordering_random_two_ints(const int32_t size, std::mt19937_64 &rng, int32_t &p, int32_t &q)
{
    std::uniform_int_distribution<> dist(0, size - 1);
    p = dist(rng);
    q = dist(rng);
    if(p > q)
    {
        int32_t temp = p;
        p = q;
        q = temp;
    }
}

void ordering_genome_swap(ORDERING_TIG *tigs, int32_t a, int32_t b)
{
    ORDERING_TIG temp = tigs[a];
    tigs[a] = tigs[b];
    tigs[b] = temp;
}

void ordering_genome_mutate(void *user, std::mt19937_64 &rng)
{
    ORDERING_GENOME *genome = reinterpret_cast<ORDERING_GENOME *>(user);
    // Mutate a Tour by applying by inversion or insertion
    std::uniform_real_distribution<double> rate(0.0, 1.0);
    double rd = rate(rng);
    if(rd < 0.2)
    {
        // MutPermute permutes two genes at random n times
        if(genome->tig_size <= 1)
        {
            return;
        }
        int32_t p, q;
        ordering_random_two_ints(static_cast<int32_t>(genome->tig_size), rng, p, q);
        if(p != q)
        {
            ordering_genome_swap(genome->tigs, p, q);
        }
    }
    else if(rd < 0.4)
    {
        // MutSplice splits a genome in 2 and glues the pieces back together in
        // reverse order.
        std::uniform_int_distribution<> dist(0, static_cast<int32_t>(genome->tig_size) - 1);
        int32_t k = dist(rng);
        int32_t tig_size = static_cast<int32_t>(genome->tig_size);
        size_t total_bytes = sizeof(ORDERING_TIG) * genome->tig_size;
        ORDERING_TIG *updated = static_cast<ORDERING_TIG *>(malloc(total_bytes));
        for(int32_t i=0; i<k; ++i)
        {
            updated[tig_size-k+i] = genome->tigs[i];
        }
        for(int32_t i=k; i<tig_size; ++i)
        {
            updated[i-k] = genome->tigs[i];
        }
        for(int32_t i=0; i<tig_size; ++i)
        {
            genome->tigs[i] = updated[i];
        }
        free(updated);
    }
    else if(rd < 0.7)
    {
        // MutInsertion applies insertion operation on the genome
        // Choose two points on the genome
        int32_t p, q;
        ordering_random_two_ints(static_cast<int32_t>(genome->tig_size), rng, p, q);
        if(p == q)
        {
            return;
        }
        std::uniform_real_distribution<double> rate(0.0, 1.0);
        if(rate(rng) < 0.5)
        {
            // Pop q and insert to p position
            ORDERING_TIG cq = genome->tigs[q];
            for(int32_t i=q; i>p; --i)
            {
                genome->tigs[i] = genome->tigs[i-1];
            }
            genome->tigs[p] = cq;
        }
        else
        {
            // Pop p and insert to q position
            ORDERING_TIG cp = genome->tigs[p];
            for(int32_t i=p; i<q; ++i)
            {
                genome->tigs[i] = genome->tigs[i+1];
            }
            genome->tigs[q] = cp;
        }
    }
    else
    {
        // MutInversion applies inversion operation on the genome
        // Choose two points on the genome
        int32_t p, q;
        ordering_random_two_ints(static_cast<int32_t>(genome->tig_size), rng, p, q);
        if(p == q)
        {
            return;
        }
        // Swap within range
        while(p < q)
        {
            ordering_genome_swap(genome->tigs, p, q);
            ++p;
            --q;
        }
    }
}

double ordering_genome_evaluate(void *user)
{
    ORDERING_GENOME *genome = reinterpret_cast<ORDERING_GENOME *>(user);
    uint64_t size = genome->tig_size;
    std::unordered_map<int32_t, double> mid;
    double cum_sum = 0.0;
    for(uint64_t i=0; i<size; ++i)
    {
        double t_size = static_cast<double>(genome->tigs[i].size);
        mid.insert(std::make_pair(genome->tigs[i].idx, cum_sum + t_size / 2.0));
        cum_sum += t_size;
    }

    double score = 0.0;
    for(auto &a_iter: *(genome->edges))
    {
        int32_t a = a_iter.first;
        double a_mid = mid.find(a)->second;
        for(auto &b_iter: a_iter.second)
        {
            int32_t b = b_iter.first, n_links = b_iter.second;
            double dist = mid.find(b)->second - a_mid;
            // This serves two purposes:
            // 1. Break earlier reduces the amount of calculation
            // 2. Ignore distant links so that telomeric regions don't come
            //    to be adjacent (based on Ler0 data)
            if(dist > LIMIT)
            {
                break;
            }
            // We are looking for maximum
            score -= double(n_links) / dist;
        }
    }
    return score;
}

void ordering_genome_crossover(void*, void *, std::mt19937_64 &)
{
    //Leave empty.
}

GA_GENOME_OP ordering_genome_op =
{
    ordering_genome_create,
    ordering_genome_free,
    ordering_genome_clone,
    ordering_genome_mutate,
    ordering_genome_evaluate,
    ordering_genome_crossover
};

ORDERING_TIG *ordering_init(ORDERING_INFO &info, int32_t *order)
{
    ORDERING_TIG *tigs = static_cast<ORDERING_TIG *>(malloc(sizeof(ORDERING_TIG) * info.contig_size));
    if(order)
    {
        time_error(-1, "Not implemented for ordered init.");
    }
    else
    {
        int32_t i = 0;
        for(auto iter: info.contigs)
        {
            tigs[i].idx = iter.first;
            tigs[i].size = iter.second.length;
            ++i;
        }
    }
    return tigs;
}

uint64_t ordering_optimize_get_int63(void *user)
{
    std::mt19937_64 *ordering_rng = reinterpret_cast<std::mt19937_64 *>(user);
    return (*ordering_rng)();
}

void ordering_callback(GA &ga)
{
    ORDERING_INFO *info = reinterpret_cast<ORDERING_INFO *>(ga.config.user);
    uint64_t gen = ga.generations;
    double current_best = -ga.hall_of_fame.best.fitness;
    if(current_best > info->best)
    {
        info->best = current_best;
        info->updated = gen;
    }
    if(gen % 500 == 0)
    {
        printf("%llu\t%lf\t%lf\n", gen, current_best, info->best);
    }
}

bool ordering_earlystop(GA &ga)
{
    ORDERING_INFO *info = reinterpret_cast<ORDERING_INFO *>(ga.config.user);
    return (ga.generations - info->updated) > ga.config.num_of_gens;
}

void ordering_optimize_phase(int32_t phase, int npop, int ngen, double mutapb,
                             ORDERING_INFO &info)
{
    // GARun set up the Genetic Algorithm and run it
    GA_CONFIG ga_config = hmr_ga_default_config();
    ga_config.num_of_pops = 1;
    ga_config.num_of_gens = 1000000;
    ga_config.pop_size = npop;
    //Configure the model.
    GA_SELECTOR_SELTOURNAMENT selector_user;
    selector_user.num_of_contestants = 3;
    GA_SELECTOR model_selector { Ga_Selector_SelTournament, &selector_user };
    GA_MODEL_GENERATIONAL model_user { model_selector, mutapb, 0.7 };
    ga_config.model.ops = Ga_Model_Generational;
    ga_config.model.model_user = &model_user;
    //Update the rng.
    std::mt19937_64 ordering_rng((std::random_device()()));
    ga_config.rng.user = &ordering_rng;
    ga_config.rng.get_int63 = ordering_optimize_get_int63;
    info.best = std::numeric_limits<double>::min();
    info.updated = 0;
    //Additional bookkeeping per generation
    ga_config.callback = ordering_callback;
    ga_config.earlystop = ordering_earlystop;
    ga_config.user = &info;
    //Start the GA.
    hmr_ga_minimize(ga_config, ordering_genome_op);
}
