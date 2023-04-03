#include <cassert>

#include "hmr_global.hpp"
#include "hmr_ui.hpp"

#include "hmr_ga.hpp"

std::vector<uint64_t> random_ints_in_range(uint64_t k, int64_t min, int64_t max, std::mt19937_64 &rng)
{
    std::vector<uint64_t> ints;
    ints.reserve(k);
    for(uint64_t i=0; i<k; ++i)
    {
        ints.push_back(i + min);
    }
    for(uint64_t i=static_cast<int64_t>(k), i_max = max-min; i<i_max; ++i)
    {
        std::uniform_int_distribution<> distr(0, i);
        uint64_t j = distr(rng);
        if(j < k)
        {
            ints[j] = i + min;
        }
    }
    return ints;
}

void sample_ints(const std::vector<uint64_t> &ints, uint64_t k, std::mt19937_64 &rng,
                 std::vector<uint64_t> &sample, std::vector<uint64_t> &idxs)
{
    if(k > ints.size())
    {
        time_error(-1, "Cannot sample %llu elements from array of length %llu", k, ints.size());
    }
    idxs = random_ints_in_range(k, 0, ints.size(), rng);
    sample.reserve(k);
    for(uint64_t i=0; i<k; ++i)
    {
        sample.push_back(ints[idxs[i]]);
    }
}

std::vector<uint64_t> generate_range_ints(uint64_t n)
{
    std::vector<uint64_t> ints;
    ints.reserve(n);
    for(uint64_t i=0; i<n; ++i)
    {
        ints.push_back(i);
    }
    return ints;
}


GA_GENOME *hmr_ga_genome_create(GA_GENOME_OP &ops, std::mt19937_64 &rng, void *user)
{
    GA_GENOME *dst = static_cast<GA_GENOME *>(malloc(sizeof(GA_GENOME)));
    dst->data = ops.create(rng, user);
    dst->ref = 1;
    return dst;
}

void hmr_ga_genome_free(GA_GENOME *genome, GA_GENOME_OP &ops)
{
    --genome->ref;
    if(genome->ref)
    {
        return;
    }
    if(genome->data)
    {
        ops.free(genome->data);
    }
    free(genome);
}

GA_GENOME *hmr_ga_genome_clone(GA_GENOME *src, GA_GENOME_OP &ops)
{
    GA_GENOME *dst = static_cast<GA_GENOME *>(malloc(sizeof(GA_GENOME)));
    dst->data = ops.clone(src->data);
    dst->ref = 1;
    return dst;
}

GA_GENOME *hmr_ga_genome_ref(GA_GENOME *genome)
{
    ++genome->ref;
    return genome;
}

double hmr_ga_genome_evaluate(GA_GENOME *genome, GA_GENOME_OP &ops)
{
    return ops.evaluate(genome->data);
}

void hmr_ga_genome_mutate(GA_GENOME *genome, GA_GENOME_OP &ops, std::mt19937_64 &rng)
{
    ops.mutate(genome->data, rng);
}

void hmr_ga_genome_crossover(GA_GENOME *src, GA_GENOME *dst, GA_GENOME_OP &ops, std::mt19937_64 &rng)
{
    ops.crossover(src->data, dst->data, rng);
}

void hmr_ga_individual_init_empty(GA_INDIVIDUAL &indi)
{
    indi.genome = NULL;
    indi.evaluated = false;
    indi.fitness = std::numeric_limits<double>::max();
    indi.id = 0;
}

void hmr_ga_individual_init(GA_INDIVIDUAL &indi, GA_GENOME *genome, std::mt19937_64 &rng)
{
    indi.genome = genome;
    indi.evaluated = false;
    indi.id = rng();
    indi.fitness = std::numeric_limits<double>::max();
}

void hmr_ga_individual_ref(GA_INDIVIDUAL &dst, GA_INDIVIDUAL &src)
{
    dst.genome = hmr_ga_genome_ref(src.genome);
    dst.evaluated = src.evaluated;
    dst.fitness = src.fitness;
    dst.id = src.id;
}

void hmr_ga_individual_clone(GA_INDIVIDUAL &dst, GA_INDIVIDUAL &src, std::mt19937_64 &rng,
                             GA_GENOME_OP &ops)
{
    dst.fitness = src.fitness;
    dst.evaluated = src.evaluated;
    dst.genome = hmr_ga_genome_clone(src.genome, ops);
    dst.id = rng();
}

void hmr_ga_individual_evaluate(GA_INDIVIDUAL &indi, GA_GENOME_OP &ops)
{
    if(indi.evaluated)
    {
        return;
    }
    indi.fitness = hmr_ga_genome_evaluate(indi.genome, ops);
    indi.evaluated = true;
}

double hmr_ga_individual_get_fitness(GA_INDIVIDUAL &indi, GA_GENOME_OP &ops)
{
    hmr_ga_individual_evaluate(indi, ops);
    return indi.fitness;
}

void hmr_ga_individual_mutate(GA_INDIVIDUAL &indi, GA_GENOME_OP &ops, std::mt19937_64 &rng)
{
    hmr_ga_genome_mutate(indi.genome, ops, rng);
    indi.evaluated = false;
}

void hmr_ga_individual_crossover(GA_INDIVIDUAL &indi, GA_INDIVIDUAL &mate,
                                 GA_GENOME_OP &ops, std::mt19937_64 &rng)
{
    hmr_ga_genome_crossover(indi.genome, mate.genome, ops, rng);
    indi.evaluated = false;
    mate.evaluated = false;
}

GA_INDIVIDUALS hmr_ga_individuals_create_empty(uint64_t n)
{
    GA_INDIVIDUALS indis;
    indis.size = n;
    indis.d = static_cast<GA_INDIVIDUAL *>(malloc(sizeof(GA_INDIVIDUAL) * n));
    assert(indis.d);
    for(uint64_t i=0; i<n; ++i)
    {
        hmr_ga_individual_init_empty(indis.d[i]);
    }
    return indis;
}

GA_INDIVIDUALS hmr_ga_individuals_create(uint64_t n, GA_GENOME_OP &ops,
                                         std::mt19937_64 &rng, void *user)
{
    GA_INDIVIDUALS indis;
    indis.size = n;
    indis.d = static_cast<GA_INDIVIDUAL *>(malloc(sizeof(GA_INDIVIDUAL) * n));
    assert(indis.d);
    for(uint64_t i=0; i<n; ++i)
    {
        hmr_ga_individual_init(indis.d[i], hmr_ga_genome_create(ops, rng, user), rng);
    }
    return indis;
}

void hmr_ga_individuals_free(GA_INDIVIDUALS &indis, GA_GENOME_OP &ops)
{
    for(uint64_t i=0; i<indis.size; ++i)
    {
        hmr_ga_genome_free(indis.d[i].genome, ops);
    }
    free(indis.d);
}

GA_INDIVIDUALS hmr_ga_individuals_clone(GA_INDIVIDUALS &src, std::mt19937_64 &rng, GA_GENOME_OP &ops)
{
    GA_INDIVIDUALS dst = hmr_ga_individuals_create_empty(src.size);
    for(size_t i=0; i<src.size; ++i)
    {
        hmr_ga_individual_clone(dst.d[i], src.d[i], rng, ops);
    }
    return dst;
}

void hmr_ga_individuals_evaluate(GA_INDIVIDUALS &indis, GA_GENOME_OP &ops)
{
    for(uint64_t i=0; i<indis.size; ++i)
    {
        hmr_ga_individual_evaluate(indis.d[i], ops);
    }
}

void hmr_ga_individuals_mutate(GA_INDIVIDUALS &indis, GA_GENOME_OP &ops,
                               const double mut_rate, std::mt19937_64 &rng)
{
    std::uniform_real_distribution<double> rate(0.0, 1.0);
    for(uint64_t i=0; i<indis.size; ++i)
    {
        if(rate(rng) < mut_rate)
        {
            hmr_ga_individual_mutate(indis.d[i], ops, rng);
        }
    }
}

void hmr_ga_individuals_sort_by_fitness(GA_INDIVIDUALS &indis)
{
    std::sort(indis.d, indis.d + indis.size, std::less<GA_INDIVIDUAL>());
}

GA_POPULATION hmr_ga_population_create(uint64_t size, GA_GENOME_OP &ops, GA_RNG &rng, void *user)
{
    GA_POPULATION pop;
    pop.rng = std::mt19937_64(rng.get_int63(rng.user));
    pop.individuals = hmr_ga_individuals_create(size, ops, pop.rng, user);
    pop.generations = 0;
    pop.id = static_cast<uint32_t>(pop.rng());
    return pop;
}

void hmr_ga_population_free(GA_POPULATION &pop, GA_GENOME_OP &ops)
{
    hmr_ga_individuals_free(pop.individuals, ops);
}

GA_INDIVIDUALS hmr_ga_seltournament_apply(uint64_t n, GA_INDIVIDUALS indis,
                                          std::mt19937_64 &rng, GA_GENOME_OP &ops, void *sel_user)
{
    GA_SELECTOR_SELTOURNAMENT *selector = reinterpret_cast<GA_SELECTOR_SELTOURNAMENT *>(sel_user);
    if((indis.size - n < selector->num_of_contestants - 1) || indis.size < n)
    {
        time_error(-1, "not enough individuals to select %llu with NContestants"
                       " = %llu, have %zu individuals and need at least %zu",
                   n, selector->num_of_contestants, indis.size,
                   selector->num_of_contestants + n - 1);
    }
    //Create empty winners, the reference will be change.
    GA_INDIVIDUALS winners = hmr_ga_individuals_create_empty(n);
    std::vector<uint64_t> not_selected_idxs = generate_range_ints(indis.size);
    for(uint64_t i=0; i<n; ++i)
    {
        std::vector<uint64_t> contestants;
        std::vector<uint64_t> idxs;
        sample_ints(not_selected_idxs, selector->num_of_contestants, rng, contestants, idxs);
        uint64_t winner_idx = idxs[0];
        for(uint64_t j=0; j<contestants.size(); ++j)
        {
            uint64_t idx = contestants[j];
            hmr_ga_individual_ref(winners.d[i], indis.d[idx]);
            winner_idx = idxs[j];
        }
        // Ban the winner from re-participating
        not_selected_idxs.erase(not_selected_idxs.begin() + winner_idx);
    }
    GA_INDIVIDUALS winner_copy = hmr_ga_individuals_clone(winners, rng, ops);
    hmr_ga_individuals_free(winners, ops);
    return winner_copy;
}


void hmr_ga_seltournament_validation(void *sel_user)
{
    GA_SELECTOR_SELTOURNAMENT *selector = reinterpret_cast<GA_SELECTOR_SELTOURNAMENT *>(sel_user);
    if(selector->num_of_contestants < 1)
    {
        time_error(-1, "NContestants should be higher than 0");
    }
}

GA_INDIVIDUALS hmr_ga_model_genearte_offsprings(uint64_t n, GA_POPULATION &pop, GA_SELECTOR &selector, GA_GENOME_OP &ops, double cross_rate)
{
    GA_INDIVIDUALS offsprings = hmr_ga_individuals_create_empty(n);
    uint64_t i=0;
    std::uniform_real_distribution<double> rate(0.0, 1.0);
    while(i < n)
    {
        //Select 2 parents.
        GA_INDIVIDUALS selected = selector.ops.apply(2, pop.individuals, pop.rng, ops, selector.selector_user);
        // Generate 2 offsprings from the parents
        if(rate(pop.rng) < cross_rate)
        {
            hmr_ga_individual_crossover(selected.d[0], selected.d[1], ops, pop.rng);
        }
        if(i < n)
        {
            hmr_ga_individual_ref(offsprings.d[i], selected.d[0]);
            ++i;
        }
        if(i < n)
        {
            hmr_ga_individual_ref(offsprings.d[i], selected.d[1]);
            ++i;
        }
        hmr_ga_individuals_free(selected, ops);
    }
    return offsprings;
}

void hmr_ga_model_generational_apply(GA_POPULATION &pop, GA_GENOME_OP &ops, void *model_user)
{
    GA_MODEL_GENERATIONAL *model = reinterpret_cast<GA_MODEL_GENERATIONAL *>(model_user);
    //Generate as many offsprings as there are of individuals in the current population.
    GA_INDIVIDUALS offsprings = hmr_ga_model_genearte_offsprings(pop.individuals.size, pop, model->selector, ops, model->cross_rate);
    if(model->mut_rate > 0)
    {
        hmr_ga_individuals_mutate(offsprings, ops, model->mut_rate, pop.rng);
    }
    //Replace the original population individuals.
    hmr_ga_individuals_free(pop.individuals, ops);
    pop.individuals = offsprings;
}

void hmr_ga_model_generational_validate(void *model_user)
{
    GA_MODEL_GENERATIONAL *model = reinterpret_cast<GA_MODEL_GENERATIONAL *>(model_user);
    // Check the selection method presence
    if(!model->selector.ops.validate)
    {
        time_error(-1, "Model Generational does not provided a selector validator.");
    }
    assert(model->selector.ops.validate);
    model->selector.ops.validate(model->selector.selector_user);
    if(model->mut_rate < 0 || model->mut_rate > 1)
    {
        time_error(-1, "Model Generational has an invalid mut rate.");
    }
    if(model->cross_rate < 0 || model->cross_rate > 1)
    {
        time_error(-1, "Model Generational has an invalid cross rate.");
    }
}

void hmr_ga_hof_init(GA_HOF &hof, size_t hof_size)
{
    hof.buffer.reserve(hof_size);
    hof.size = hof_size;
    hof.queue = GA_HOF_QUEUE(std::greater<GA_INDIVIDUAL>(), std::move(hof.buffer));
    hmr_ga_individual_init_empty(hof.best);
}

void hmr_ga_hof_update(GA_HOF &hof, GA_INDIVIDUALS &indis, std::mt19937_64 &rng, GA_GENOME_OP &ops)
{
    // Start by finding the current best Individual.
    for(size_t i=0, i_max = hMin(indis.size, hof.size); i<i_max; ++i)
    {
        GA_INDIVIDUAL &indi = indis.d[i];
        // Find if and where the Individual should fit in the hall of fame
        if(hof.queue.empty() || indi.fitness < hof.queue.top().fitness)
        {
            //Update the best record first.
            GA_INDIVIDUAL fame_indi;
            hmr_ga_individual_clone(fame_indi, indi, rng, ops);
            //Check the best fit.
            if(fame_indi.fitness < hof.best.fitness)
            {
                hof.best = fame_indi;
            }
            if(hof.queue.size() == hof.size)
            {
                hmr_ga_genome_free(hof.queue.top().genome, ops);
                hof.queue.pop();
            }
            hof.queue.push(fame_indi);
        }
    }
}

void hmr_ga_init(GA &ga, GA_GENOME_OP &ops)
{
    GA_CONFIG &config = ga.config;
    //Check config validation.
    if(config.num_of_pops == 0)
    {
        time_error(-1, "NPops has to be strictly higher than 0");
    }
    if(config.pop_size == 0)
    {
        time_error(-1, "PopSize has to be strictly higher than 0");
    }
    if(config.num_of_gens == 0)
    {
        time_error(-1, "NGenerations has to be strictly higher than 0");
    }
    if(config.hof_size == 0)
    {
        time_error(-1, "HofSize has to be strictly higher than 0");
    }
    if(!config.model.ops.validate)
    {
        time_error(-1, "model has to be provided");
    }
    assert(config.model.ops.validate);
    config.model.ops.validate(config.model.model_user);
    // Reset counters
    ga.generations = 0;
    // Create the initial Populations
    ga.populations.size = config.num_of_pops;
    ga.populations.d = static_cast<GA_POPULATION *>(malloc(sizeof(GA_POPULATION) * ga.populations.size));
    for(uint64_t i=0; i<ga.populations.size; ++i)
    {
        ga.populations.d[i] = hmr_ga_population_create(config.pop_size, ops, config.rng, config.user);
        hmr_ga_individuals_evaluate(ga.populations.d[i].individuals, ops);
        hmr_ga_individuals_sort_by_fitness(ga.populations.d[i].individuals);
    }
    //Initiallize the hall of fame.
    hmr_ga_hof_init(ga.hall_of_fame, config.hof_size);
    for(uint64_t i=0; i<ga.populations.size; ++i)
    {
        GA_POPULATION &pop = ga.populations.d[i];
        hmr_ga_hof_update(ga.hall_of_fame, pop.individuals, pop.rng, ops);
    }
    // Execute the callback if it has been set.
    if(config.callback)
    {
        config.callback(ga);
    }
}

void hmr_ga_evolve(GA &ga, GA_GENOME_OP &ops)
{
    ++ga.generations;
    GA_CONFIG &config = ga.config;

    for(uint64_t i=0; i<ga.populations.size; ++i)
    {
        //Apply the evolution model to the entire population
        GA_POPULATION &pop = ga.populations.d[i];
        config.model.ops.apply(pop, ops, config.model.model_user);
        //Evaluate and sort.
        hmr_ga_individuals_evaluate(pop.individuals, ops);
        hmr_ga_individuals_sort_by_fitness(pop.individuals);
        // Record time spent evolving.
        ++pop.generations;
    }
    // Update HallOfFame
    for(uint64_t i=0; i<ga.populations.size; ++i)
    {
        GA_POPULATION &pop = ga.populations.d[i];
        hmr_ga_hof_update(ga.hall_of_fame, pop.individuals, pop.rng, ops);
    }
    // Execute the callback if it has been set.
    if(config.callback)
    {
        config.callback(ga);
    }
}

GA_CONFIG hmr_ga_default_config()
{
    GA_CONFIG config;
    config.num_of_pops = 1;
    config.pop_size = 30;
    config.num_of_gens = 50;
    config.hof_size = 1;
    config.model = GA_MODEL { GA_MODEL_OPS {NULL, NULL}, NULL };

    config.parallel_init = false;
    config.parallel_eval = false;
    config.callback = NULL;
    config.earlystop = NULL;
    config.rng = GA_RNG { NULL, NULL },
    config.user = NULL;
    return config;
}


void hmr_ga_minimize(const GA_CONFIG &config, GA_GENOME_OP &ops)
{
    GA ga;
    ga.config = config;
    time_print("Initialize GA parameters...");
    hmr_ga_init(ga, ops);
    // Go through the generations
    for(uint64_t i=0; i<ga.config.num_of_gens; ++i)
    {
        //Check for early stopping.
        if(ga.config.earlystop && ga.config.earlystop(ga))
        {
            return;
        }
        hmr_ga_evolve(ga, ops);
    }
}
