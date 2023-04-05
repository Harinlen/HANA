#include <cassert>
#include <numeric>
#include <limits>
#include <random>
#include <cmath>

#include "hmr_global.hpp"
#include "hmr_parallel.hpp"
#include "hmr_ui.hpp"

#include "ordering_descent.hpp"

constexpr double LIMIT = 10000000;
const double LimitLog = log(LIMIT);

typedef struct ORDERING_TIG
{
    int32_t index;
    int32_t length;
} ORDERING_TIG;

typedef struct ORDERING_EVA
{
    ORDERING_TIG *seq;
    double score;
    bool evaluated;
} ORDERING_EVA;

ORDERING_TIG *ordering_init(ORDERING_INFO &info, bool shuffle, const int32_t *order)
{
    ORDERING_TIG *tigs = static_cast<ORDERING_TIG *>(malloc(sizeof(ORDERING_TIG) * info.contig_size));
    if(order)
    {
        for(int32_t i=0; i<info.contig_size; ++i)
        {
            tigs[i].index = order[i];
            auto iter = info.contigs.find(tigs[i].index);
            if(iter == info.contigs.end())
            {
                time_error(-1, "Failed to find contig %d", tigs[i].index);
            }
            tigs[i].length = iter->second.length;
        }
    }
    else
    {
        int32_t i = 0;
        for(auto &iter: info.contigs)
        {
            tigs[i].index = iter.first;
            tigs[i].length = iter.second.length;
            ++i;
        }
        if(shuffle)
        {
            uint64_t seed = std::mt19937_64((std::random_device()()))();
            std::shuffle(tigs, tigs+info.contig_size,
                         std::default_random_engine(static_cast<unsigned>(seed)));
        }
    }
    return tigs;
}

typedef struct ORDERING_EA
{
    int32_t num_of_seqs;
    double mut_rate;
    uint64_t generations;
    std::mt19937_64 rng;

    ORDERING_TIG *buffer1, *buffer2; // Sequence buffer pool.
    ORDERING_EVA *eva1, *eva2;

    ORDERING_TIG *best_seq;
    double best_score;
    uint64_t best_gen;

    bool reversed;
    ORDERING_EVA *current_eva, *target_eva;
    ORDERING_TIG *target_buffer;
} ORDERING_EA;

int32_t ordering_ea_get_links(int32_t a, int32_t b, const ORDERING_COUNTS &edges)
{
    if(a > b)
    {
        int32_t temp = a;
        a = b;
        b = temp;
    }
    const auto &a_iter = edges.find(a);
    if(a_iter == edges.end())
    {
        return 0;
    }
    const auto &b_iter = a_iter->second.find(b);
    return b_iter == a_iter->second.end() ? 0 : b_iter->second;
}

double ordering_evaluate_sequence(ORDERING_TIG *seq, int32_t seq_length, const ORDERING_COUNTS &edges)
{
    double *mid = static_cast<double *>(malloc(sizeof(double) * seq_length));
    double cum_sum = 0.0;
    for(int32_t i=0; i<seq_length; ++i)
    {
        double t_size = static_cast<double>(seq[i].length);
        mid[i] = cum_sum + t_size / 2.0;
        cum_sum += t_size;
    }
    //Calculate the score of the sequence.
    double score = 0.0;
    for(int32_t i=0; i<seq_length - 1; ++i)
    {
        int32_t a = seq[i].index;
        double a_mid = mid[i];
        for(int32_t j=i+1; j<seq_length; ++j)
        {
            int32_t b = seq[j].index,
                    n_links = ordering_ea_get_links(a, b, edges);
            double dist = mid[j] - a_mid;
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
    free(mid);
    return score;
}

void ordering_evaluate_calc(ORDERING_EVA *eva, int32_t seqs_count, int32_t seq_length,
                            const ORDERING_COUNTS &edges)
{
    //Loop and check all the evaluation state.
    ;
    for(int32_t i=0; i<seqs_count; ++i)
    {
        if(!eva[i].evaluated)
        {
            eva[i].score = ordering_evaluate_sequence(eva[i].seq, seq_length, edges);
            eva[i].evaluated = true;
        }
    }
}

void ordering_evaluate_sort(ORDERING_EVA *eva, int32_t seqs_count)
{
    std::sort(eva, eva + seqs_count,
              [](const ORDERING_EVA &lhs, const ORDERING_EVA &rhs)
    {
        return lhs.score > rhs.score;
    });
}

typedef std::vector<int32_t> ORDERING_IDS;
ORDERING_IDS ordering_select_contestants(int32_t n, ORDERING_IDS &candidates, std::mt19937_64 &rng)
{
    std::unordered_set<int32_t> chosen;
    std::uniform_int_distribution<> dist(0, candidates.size() - 1);
    for(int32_t i=0; i<n; ++i)
    {
        int32_t j = dist(rng);
        while (chosen.count(j) > 0)
        {
            j = dist(rng);
        }
        chosen.insert(j);
    }
    return ORDERING_IDS(chosen.begin(), chosen.end());
}

// SelTournament samples individuals through tournament selection. The
// tournament is composed of randomly chosen individuals. The winner of the
// tournament is the chosen individual with the lowest fitness. The obtained
// individuals are all distinct, in other words there are no repetitions.
ORDERING_IDS ordering_select_n_parents(int32_t n, int32_t num_of_contestants,
                                       ORDERING_EVA *parents, int32_t seqs_count, std::mt19937_64 &rng)
{
    // Check that the number of individuals is large enough
    if(seqs_count - n < num_of_contestants - 1 || seqs_count < n)
    {
        time_error(-1, "Not enough contestants.");
    }
    ORDERING_IDS winners;
    ORDERING_IDS candidates(seqs_count);
    winners.reserve(n);
    std::iota(candidates.begin(), candidates.end(), 0);
    for(int32_t i=0; i<n; ++i)
    {
        // Select contestants.
        std::vector<int32_t> contestants_id = ordering_select_contestants(num_of_contestants, candidates, rng);
        // Find the winner from contestants.
        int32_t winner_id = contestants_id[0];
        double winner_score = parents[winner_id].score;
        for(int32_t j=1; j<num_of_contestants; ++j)
        {
            int32_t current_id = contestants_id[j];
            double current_score = parents[current_id].score;
            if(current_score < winner_score)
            {
                winner_score = current_score;
                winner_id = current_id;
            }
        }
        //Record the winner id.
        winners.push_back(winner_id);
        // Ban the winner from re-participating.
        candidates.erase(candidates.begin() + winner_id);
    }
    return winners;
}

void ordering_clone_seq(ORDERING_EVA &dst, const ORDERING_EVA &src, size_t seq_bytes)
{
    //Copy the sequence first.
    memcpy(dst.seq, src.seq, seq_bytes);
    dst.score = src.score;
    dst.evaluated = src.evaluated;
}

void ordering_generate_offsprings(ORDERING_EVA *parents, ORDERING_EVA *offsprings,
                                  int32_t seqs_count, size_t seq_bytes, std::mt19937_64 &rng)
{
    //Loop until the offsprings are filled.
    int32_t offspring_index = 0;
    while(offspring_index < seqs_count)
    {
        //Select 2 parents from 3 contestants.
        ORDERING_IDS selected = ordering_select_n_parents(2, 3, parents, seqs_count, rng);
        //Create the offsprings.
        if(offspring_index < seqs_count)
        {
            ordering_clone_seq(offsprings[offspring_index], parents[selected[0]], seq_bytes);
            ++offspring_index;
        }
        if(offspring_index < seqs_count)
        {
            ordering_clone_seq(offsprings[offspring_index], parents[selected[1]], seq_bytes);
            ++offspring_index;
        }
    }
}

void ordering_ea_tig_random_range(std::uniform_int_distribution<> &dist,
                                   std::mt19937_64 &rng, int32_t &p, int32_t &q)
{
    p = dist(rng);
    q = dist(rng);
    if(p > q)
    {
        int32_t temp = p;
        p = q;
        q = temp;
    }
}

void ordering_swap(ORDERING_TIG *seq, int32_t p, int32_t q)
{
    ORDERING_TIG temp = seq[p];
    seq[p] = seq[q];
    seq[q]= temp;
}

void ordering_mutate(ORDERING_EVA &candidate, int32_t seq_length, std::mt19937_64 &rng)
{
    //Reset the evaluate state.
    candidate.evaluated = false;
    // Mutate a Tour by applying by inversion or insertion
    std::uniform_real_distribution<> rate(0.0, 1.0);
    std::uniform_int_distribution<> index_rng(0, seq_length - 1);
    double mutate_rate = rate(rng);
    if(mutate_rate < 0.2)
    {
        //MutPermute(r, rng)
        if(seq_length < 2)
        {
            return;
        }
        // Choose two points on the genome, swap them.
        int32_t p, q;
        ordering_ea_tig_random_range(index_rng, rng, p, q);
        if(p == q)
        {
            return;
        }
        ordering_swap(candidate.seq, p, q);
    }
    else if(mutate_rate < 0.4)
    {
        //MutSplice(r, rng)
        int32_t k = index_rng(rng);
        //Construct [k:][:k].
        std::vector<ORDERING_TIG> buffer;
        buffer.reserve(seq_length);
        for(int32_t i=k; i<seq_length; ++i)
        {
            buffer.push_back(candidate.seq[i]);
        }
        for(int32_t i=0; i<k; ++i)
        {
            buffer.push_back(candidate.seq[i]);
        }
        for(int32_t i=0; i<seq_length; ++i)
        {
            candidate.seq[i] = buffer[i];
        }
    }
    else if(mutate_rate < 0.7)
    {
        //MutInsertion(r, rng)
        // Choose two points on the genome
        int32_t p, q;
        ordering_ea_tig_random_range(index_rng, rng, p, q);
        //Decide the direction.
        if(rate(rng) < 0.5)
        {
            // Pop q and insert to p position
            auto cq = candidate.seq[q];
            for(int32_t i=q; i>p; --i)
            {
                candidate.seq[i] = candidate.seq[i-1];
            }
            candidate.seq[p] = cq;
        }
        else
        {
            auto cp = candidate.seq[p];
            // Move cq to the back and push everyone left
            for(int32_t i=p; i<q; ++i)
            {
                candidate.seq[i] = candidate.seq[i+1];
            }
            candidate.seq[q] = cp;
        }
    }
    else
    {
        //MutInversion(r, rng)
        // MutInversion applies inversion operation on the genome
        // Choose two points on the genome
        int32_t p, q;
        ordering_ea_tig_random_range(index_rng, rng, p, q);
        // Swap within range
        while(p < q)
        {
            ordering_swap(candidate.seq, p, q);
            ++p;
            --q;
        }
    }
}

void ordering_update_best(ORDERING_EA &ea, size_t seq_bytes)
{
    ORDERING_EVA &best_eva = ea.current_eva[0];
    double best_score = -best_eva.score;
    if(best_score > ea.best_score)
    {
        ea.best_score = best_score;
        ea.best_gen = ea.generations;
        memcpy(ea.best_seq, best_eva.seq, seq_bytes);
    }
}

void ordering_optimize_phase(int32_t phase, int npop, int ngen, double mutapb,
                             ORDERING_INFO &info)
{
    ORDERING_EA ea;
    ea.num_of_seqs = npop;
    ea.mut_rate = mutapb;
    ea.generations = 0;
    ea.rng = std::mt19937_64(std::random_device()());
    //Allocate memory for the buffer.
    size_t seq_bytes = sizeof(ORDERING_TIG) * info.contig_size,
            buffer_bytes = ea.num_of_seqs * seq_bytes;
    ea.buffer1 = static_cast<ORDERING_TIG *>(malloc(buffer_bytes));
    ea.buffer2 = static_cast<ORDERING_TIG *>(malloc(buffer_bytes));
    size_t eva_bytes = sizeof(ORDERING_EVA) * ea.num_of_seqs;
    ea.eva1 = static_cast<ORDERING_EVA *>(malloc(eva_bytes));
    ea.eva2 = static_cast<ORDERING_EVA *>(malloc(eva_bytes));
    //Activate sequence 1.
    ea.reversed = false;
    ea.current_eva = ea.eva1;
    ea.target_eva = ea.eva2;
    ea.target_buffer = ea.buffer2;
    std::uniform_real_distribution<double> rate(0.0, 1.0);
    //Assign the sequence pointer.
    for(int32_t i=0; i<ea.num_of_seqs; ++i)
    {
        //Set and initialize the buffer1.
        off_t seq_offset = i * info.contig_size;
        ea.eva1[i].seq = ea.buffer1 + seq_offset;
        //Initialize the sequences.
        memcpy(ea.eva1[i].seq, info.init_genome, seq_bytes);
        ea.eva1[i].evaluated = false;
        ea.eva1[i].score = std::numeric_limits<double>::max();
    }
    //Evaluate the current buffer.
    ordering_evaluate_calc(ea.current_eva, ea.num_of_seqs, info.contig_size, info.edges);
    ordering_evaluate_sort(ea.current_eva, ea.num_of_seqs);
    //Initialize the best records.
    ea.best_score = -std::numeric_limits<double>::max();
    ea.best_gen = 0;
    ea.best_seq = static_cast<ORDERING_TIG *>(malloc(seq_bytes));
    //Update the best result.
    ordering_update_best(ea, seq_bytes);
    //Run EA algorithm with the configuration above.
    for(int32_t g=0; g<1000000; ++g)
    {
        //Convergence criteria.
        if(ea.generations - ea.best_gen > static_cast<uint64_t>(ngen))
        {
            //Meet the requirement, so we quit.
            break;
        }
        //Increase the generations.
        ++ea.generations;
        std::mt19937_64 gen_rng(ea.rng());
        //Reset the target eva positions.
        for(int32_t i=0; i<ea.num_of_seqs; ++i)
        {
            ea.target_eva[i].seq = ea.target_buffer + i * info.contig_size;
        }
        //Generate the offspring at target eva matrix.
        ordering_generate_offsprings(ea.current_eva, ea.target_eva, ea.num_of_seqs, seq_bytes, gen_rng);
        if(ea.mut_rate > 0)
        {
            for(int32_t i=0; i<ea.num_of_seqs; ++i)
            {
                if(rate(gen_rng) < ea.mut_rate)
                {
                    ordering_mutate(ea.target_eva[i], info.contig_size, gen_rng);
                }
            }
        }
        //Evaluate the target buffer.
        ordering_evaluate_calc(ea.target_eva, ea.num_of_seqs, info.contig_size, info.edges);
        ordering_evaluate_sort(ea.target_eva, ea.num_of_seqs);
        //Update the best result.
        ordering_update_best(ea, seq_bytes);
        printf("g=%d\t%lf\t%lf\t%lf\t%lf\n", g, ea.target_eva[0].score, ea.target_eva[1].score, ea.target_eva[2].score, ea.best_score);
        //Swap the target and current.
        if(ea.reversed)
        {
            ea.reversed = false;
            ea.current_eva = ea.eva1;
            ea.target_eva = ea.eva2;
            ea.target_buffer = ea.buffer2;
        }
        else
        {
            ea.reversed = true;
            ea.current_eva = ea.eva2;
            ea.target_eva = ea.eva1;
            ea.target_buffer = ea.buffer1;
        }
    }
    for(int32_t i=0; i<info.contig_size; ++i)
    {
        printf("%d\t%s\n", ea.best_seq[i].index, info.contigs.find(ea.best_seq[i].index)->second.name);
    }
}
