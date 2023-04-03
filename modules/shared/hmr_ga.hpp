#ifndef HMR_GA_H
#define HMR_GA_H

#include <cstdint>
#include <queue>
#include <random>

typedef struct GA_GENOME_OP
{
    void *(*create)(std::mt19937_64 &, void *);
    void (*free)(void*);
    void *(*clone)(void*);
    void (*mutate)(void*, std::mt19937_64 &);
    double (*evaluate)(void*);
    void (*crossover)(void*, void *, std::mt19937_64 &);
} GA_GENOME_OP;

typedef struct GA_GENOME
{
    void *data;
    uint32_t ref;
} GA_GENOME;

typedef struct GA_INDIVIDUAL
{
    GA_GENOME *genome;
    double fitness;
    bool evaluated;
    uint64_t id;
    bool operator> (const GA_INDIVIDUAL &other) const
    {
        return fitness > other.fitness;
    }

    bool operator< (const GA_INDIVIDUAL &other) const
    {
        return fitness < other.fitness;
    }
} GA_INDIVIDUAL;

typedef struct GA_INDIVIDUALS
{
    GA_INDIVIDUAL *d;
    size_t size;
} GA_INDIVIDUALS;

typedef struct GA_SELECTOR_OPS
{
    GA_INDIVIDUALS (*apply)(uint64_t n, GA_INDIVIDUALS indis, std::mt19937_64 &rng, GA_GENOME_OP &ops, void *);
    void (*validate)(void *);
} GA_SELECTOR_OPS;

GA_INDIVIDUALS hmr_ga_seltournament_apply(uint64_t, GA_INDIVIDUALS, std::mt19937_64 &, GA_GENOME_OP &, void *);
void hmr_ga_seltournament_validation(void *);

typedef struct GA_SELECTOR_SELTOURNAMENT
{
    uint64_t num_of_contestants;
} GA_SELECTOR_SELTOURNAMENT;

static GA_SELECTOR_OPS Ga_Selector_SelTournament
{
    hmr_ga_seltournament_apply,
    hmr_ga_seltournament_validation
};

typedef struct GA_SELECTOR
{
    GA_SELECTOR_OPS ops;
    void *selector_user;
} GA_SELECTOR;

typedef struct GA_POPULATION
{
    GA_INDIVIDUALS individuals;
    uint64_t generations;
    uint32_t id;
    std::mt19937_64 rng;
} GA_POPULATION;

typedef struct GA_POPULATIONS
{
    GA_POPULATION *d;
    size_t size;
} GA_POPULATIONS;

typedef struct GA_MODEL_OPS
{
    void (*apply)(GA_POPULATION &, GA_GENOME_OP &, void *);
    void (*validate)(void *);
} GA_MODEL_OPS;

void hmr_ga_model_generational_apply(GA_POPULATION &, GA_GENOME_OP &, void *);
void hmr_ga_model_generational_validate(void *);

typedef struct GA_MODEL_GENERATIONAL
{
    GA_SELECTOR selector;
    double mut_rate;
    double cross_rate;
} GA_MODEL_GENERATIONAL;

static GA_MODEL_OPS Ga_Model_Generational =
{
    hmr_ga_model_generational_apply,
    hmr_ga_model_generational_validate
};

typedef struct GA_MODEL
{
    GA_MODEL_OPS ops;
    void *model_user;
} GA_MODEL;

typedef struct GA_RNG
{
    uint64_t (*get_int63)(void *);
    void *user;
} GA_RNG;

typedef std::priority_queue<GA_INDIVIDUAL, std::vector<GA_INDIVIDUAL>, std::greater<GA_INDIVIDUAL> > GA_HOF_QUEUE;

typedef struct GA_HOF
{
    std::vector<GA_INDIVIDUAL> buffer;
    GA_HOF_QUEUE queue;
    GA_INDIVIDUAL best;
    size_t size;
} GA_HOF;

typedef struct GA GA;

typedef void (*GA_CALLBACK)(GA &);
typedef bool (*GA_EARLYSTOP)(GA &);

typedef struct GA_CONFIG
{
    uint64_t num_of_pops;
    uint64_t pop_size;
    uint64_t num_of_gens;
    size_t hof_size;
    GA_MODEL model;

    bool parallel_init;
    bool parallel_eval;
    GA_CALLBACK callback;
    GA_EARLYSTOP earlystop;
    GA_RNG rng;
    void *user;

} GA_CONFIG;

typedef struct GA
{
    GA_CONFIG config;
    GA_POPULATIONS populations;
    uint64_t generations;
    GA_HOF hall_of_fame;
} GA;

GA_CONFIG hmr_ga_default_config();
void hmr_ga_minimize(const GA_CONFIG &config, GA_GENOME_OP &ops);

#endif // HMR_GA_H
