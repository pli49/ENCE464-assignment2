#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
 * poisson.c
 * Implementation of a Poisson solver with Neumann boundary conditions.
 *
 * This template handles the basic program launch, argument parsing, and memory
 * allocation required to implement the solver *at its most basic level*. You
 * will likely need to allocate more memory, add threading support, account for
 * cache locality, etc...
 *
 * BUILDING:
 * g++ -pg -o poisson poisson.c -lpthread
 *
 * [note: linking pthread isn't strictly needed until you add your
 *        multithreading code]
 *
 * TODO:
 * 1 - Read through this example, understand what it does and what it gives you
 *     to work with.
 * 2 - Implement the basic algorithm and get a correct output.
 * 3 - Add a timer to track how long your execution takes.
 * 4 - Profile your solution and identify weaknesses.
 * 5 - Improve it!
 * 6 - Remember that this is now *your* code and *you* should modify it however
 *     needed to solve the assignment.
 *
 * See the lab notes for a guide on profiling and an introduction to
 * multithreading (see also threads.c which is reference by the lab notes).
 */


// Global flag
// set to true when operating in debug mode to enable verbose logging
static bool debug = false;


/**
 * @brief Solve Poissons equation for a given cube with Neumann boundary
 * conditions on all sides.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing function.
 * @param iterations    Number of iterations to solve with.
 * @param threads       Number of threads to use for solving.
 * @param delta         Grid spacing.
 * @return double*      Solution to Poissons equation. Caller must free().
 */


typedef struct {
    int xstart;
    int xend;
    int ystart;
    int yend;
    int zstart;
    int zend;

    double* curr;
    int n;
    double *source;
    float delta;
    double* next;

} WorkerArgs;

static inline double xadd_both(double* curr, int x, int y, int z, int n) {
return curr[(z*n + y) * n + (x - 1)] + //added left 
         curr[(z*n + y) * n + (x + 1)]; //added right
}

static inline double xzero(double* curr, int y, int z, int n) {
    return 2 * curr[(z*n + y) * n + 1]; //add twice right  
}

static inline double xn_one(double* curr, int y, int z, int n) {
    return 2 * curr[(z*n + y) * n + (n - 2)];
}

static inline double yadd_both(double* curr, int x, int y, int z, int n) {
    return curr[(z*n + (y - 1)) * n + x]
            + curr[(z*n + (y + 1)) * n + x];
}

static inline double yzero(double* curr, int x, int z, int n) {
    return 2 * curr[(z*n + 1) * n + x];
}

static inline double yn_one(double* curr, int x, int z, int n) {
    return 2 * curr[(z*n + (n - 2)) * n + x];
}

static inline double zadd_both(double* curr, int x, int y, int z, int n) {
    return curr[((z - 1)*n + y) * n + x] +
            curr[((z + 1)*n + y) * n + x];
}

static inline double zzero(double* curr, int x, int y, int n) {
    return 2 * curr[(n + y) * n + x];
}

static inline double zn_one(double* curr, int x, int y, int n) {
    return 2 * curr[((n - 2)*n + y) * n + x];
}


void equation(void* pargs)
{
    WorkerArgs* args = (WorkerArgs*)pargs;

    double sum;
    int zloop_start = args->zstart;
    int yloop_start;
    int xloop_start;

    // To completely unroll these loops, since 6 boundaries, 64 combinations
    // Could do it but nah. Reckon we'd need to write code to write the code to make it worth doing that
    // Possible pro gamer move if we did though: make the threads at the beginning of the program,
    // 6 layers of if statements to determine which code it should run, then rendezvous between them
    // last one there switches current and next (using pointers to pointers). wouldn't need to create
    // and collect the threads each iteration either. Outermost loop though so probs wouldn't make a big difference

    // Either way, for medium to large n, a common case is no border in x or y,
    // check this first so that there are fewer coomparisons
    if ( args->ystart != 0 && args->yend != args->n - 1
            && args->xstart != 0 && args->xend != args->n - 1 ) {

        // z lower boundary?
        if (args->zstart == 0) {
            for (int y = args->ystart; y <= args->yend; y++) {
                for (int x = args->xstart; x <= args->xend - 1; x++) {
                    sum = xadd_both(args->curr, x, y, 0, args->n);
                    sum += yadd_both(args->curr, x, y, 0, args->n);
                    sum += zzero(args->curr, x, y, args->n);

                    sum -= args->delta * args->delta * args->source[y* args->n + x];
                    args->next[y*args->n + x] = sum / 6;
                }
            }
            zloop_start += 1;
        }

        // inner z
        for (int z = zloop_start; z < args->zend; z++) {
            for (int y = args->ystart; y <= args->yend; y++) {
                for (int x = args->xstart; x <= args->xend - 1; x++) {
                    sum = xadd_both(args->curr, x, y, 0, args->n);
                    sum += yadd_both(args->curr, x, y, 0, args->n);
                    sum += zadd_both(args->curr, x, y, z, args->n);

                    sum -= args->delta * args->delta * args->source[(z*args->n + y) * args->n + x];
                    args->next[(z*args->n + y)*args->n + x] = sum / 6;
                }
            }
        }

        // z upper boundary
        if (args->zend == args->n - 1) {
            for (int y = args->ystart; y <= args->yend; y++) {
                for (int x = args->xstart; x <= args->xend - 1; x++) {
                    sum = xadd_both(args->curr, x, y, 0, args->n);
                    sum += yadd_both(args->curr, x, y, 0, args->n);
                    sum += zn_one(args->curr, x, y, args->n);

                    sum -= args->delta * args->delta * args->source[(args->zend*args->n + y) * args->n + x];
                    args->next[(args->zend*args->n + y)*args->n + x] = sum / 6;
                }
            }
        }

        // no z upper boundary
        else {
            for (int y = args->ystart; y <= args->yend; y++) {
                for (int x = args->xstart; x <= args->xend - 1; x++) {
                    sum = xadd_both(args->curr, x, y, 0, args->n);
                    sum += yadd_both(args->curr, x, y, 0, args->n);
                    sum += zadd_both(args->curr, x, y, args->zend, args->n);

                    sum -= args->delta * args->delta * args->source[(args->zend*args->n + y) * args->n + x];
                    args->next[(args->zend*args->n + y)*args->n + x] = sum / 6;
                }
            }
        }
    }


    // at least one of x and y have a boundary. taking x and y if statements out of inner most loop
    else {
        yloop_start = (args->ystart == 0) ? 1 : args->ystart;
        xloop_start = (args->xstart == 0) ? 1 : args->xstart;

        for (int z = args->zstart; z <= args->zend; z++) {

            // if y lower boundary
            if (args->ystart == 0) {

                // if x lower boundary
                if (args->xstart == 0) {

                    sum = xzero(args->curr, 0, z, args->n);
                    sum += yzero(args->curr, 0, z, args->n);

                    // common case fist
                    if ( z == 0) {
                        sum += zzero(args->curr, 0, 0, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, 0, 0, z, args->n);
                    } else {
                        sum += zn_one(args->curr, 0, 0, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n) * args->n];
                    args->next[(z*args->n)*args->n] = sum / 6;
                }

                // x inner loop
                for (int x = xloop_start; x < args->xend; x++) {
                    sum = xadd_both(args->curr, x, 0, z, args->n);
                    sum += yzero(args->curr, x, z, args->n);

                    if ( z == 0) {
                        sum += zzero(args->curr, x, 0, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, x, 0, z, args->n);
                    } else {
                        sum += zn_one(args->curr, x, 0, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n) * args->n + x];
                    args->next[(z*args->n)*args->n + x] = sum / 6;
                }

                // if x upper boundary
                if (args->xend == args->n - 1) {
                    sum = xn_one(args->curr, 0, z, args->n);
                    sum += yzero(args->curr, args->xend, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, args->xend, 0, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, args->xend, 0, z, args->n);
                    } else {
                        sum += zn_one(args->curr, args->xend, 0, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n) * args->n + args->xend];
                    args->next[(z*args->n)*args->n + args->xend] = sum / 6;
                } else {
                    sum = xadd_both(args->curr, args->xend, 0, z, args->n);
                    sum += yzero(args->curr, args->xend, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, args->xend, 0, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, args->xend, 0, z, args->n);
                    } else {
                        sum += zn_one(args->curr, args->xend, 0, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n) * args->n + args->xend];
                    args->next[(z*args->n)*args->n + args->xend] = sum / 6;
                }
            }


            //y inner loop
            for (int y = yloop_start; y < args->yend; y++) {
                // if x lower boundary
                if (args->xstart == 0) {

                    sum = xzero(args->curr, y, z, args->n);
                    sum += yadd_both(args->curr, 0, y, z, args->n);

                    if ( z == 0) {
                        sum += zzero(args->curr, 0, y, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, 0, y, z, args->n);
                    } else {
                        sum += zn_one(args->curr, 0, y, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + y) * args->n];
                    args->next[(z*args->n + y)*args->n] = sum / 6;
                }

                // x inner loop
                for (int x = xloop_start; x < args->xend; x++) {
                    sum = xadd_both(args->curr, x, y, z, args->n);
                    sum += yadd_both(args->curr, x, y, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, x, y, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, x, y, z, args->n);
                    } else {
                        sum += zn_one(args->curr, x, y, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + y) * args->n + x];
                    args->next[(z*args->n + y)*args->n + x] = sum / 6;
                }
                // if x upper boundary
                if (args->xend == args->n - 1) {
                    sum = xn_one(args->curr, y, z, args->n);
                    sum += yadd_both(args->curr, args->xend, y, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, args->xend, y, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, args->xend, y, z, args->n);
                    } else {
                        sum += zn_one(args->curr, args->xend, y, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + y) * args->n + args->xend];
                    args->next[(z*args->n + y)*args->n + args->xend] = sum / 6;
                } else {
                    sum = xadd_both(args->curr, args->xend, y, z, args->n);
                    sum += yadd_both(args->curr, args->xend, y, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, args->xend, y, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, args->xend, y, z, args->n);
                    } else {
                        sum += zn_one(args->curr, args->xend, y, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + y) * args->n + args->xend];
                    args->next[(z*args->n + y)*args->n + args->xend] = sum / 6;
                }
            }


            // if y upper boundary

            if (args->yend == args->n - 1) {
                // if x lower boundary
                if (args->xstart == 0) {

                    sum = xzero(args->curr, args->yend, z, args->n);
                    sum += yn_one(args->curr, 0, z, args->n);

                    if ( z == 0) {
                        sum += zzero(args->curr, 0, args->yend, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, 0, args->yend, z, args->n);
                    } else {
                        sum += zn_one(args->curr, 0, args->yend, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + args->yend) * args->n];
                    args->next[(z*args->n + args->yend)*args->n] = sum / 6;
                }

                // x inner loop
                for (int x = xloop_start; x < args->xend; x++) {
                    sum = xadd_both(args->curr, x, args->yend, z, args->n);
                    sum += yn_one(args->curr, x, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, x, args->yend, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, x, args->yend, z, args->n);
                    } else {
                        sum += zn_one(args->curr, x, args->yend, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + args->yend) * args->n + x];
                    args->next[(z*args->n + args->yend)*args->n + x] = sum / 6;
                }
                // if x upper boundary
                if (args->xend == args->n - 1) {
                    sum = xn_one(args->curr, args->yend, z, args->n);
                    sum += yn_one(args->curr, args->xend, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, args->xend, args->yend, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, args->xend, args->yend, z, args->n);
                    } else {
                        sum += zn_one(args->curr, args->xend, args->yend, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + args->yend) * args->n + args->xend];
                    args->next[(z*args->n + args->yend)*args->n + args->xend] = sum / 6;
                } else {
                    sum = xadd_both(args->curr, args->xend, args->yend, z, args->n);
                    sum += yn_one(args->curr, args->xend, z, args->n);
                    
                    if ( z == 0) {
                        sum += zzero(args->curr, args->xend, args->yend, args->n);
                    } else if (z != args->n -1) {
                        sum += zadd_both(args->curr, args->xend, args->yend, z, args->n);
                    } else {
                        sum += zn_one(args->curr, args->xend, args->yend, args->n);
                    }

                    sum -= args->delta * args->delta * args->source[(z*args->n + args->yend) * args->n + args->xend];
                    args->next[(z*args->n + args->yend)*args->n + args->xend] = sum / 6;
                }
            }
        }
    }
}



double* poisson_neumann(int n, double *source, int iterations, int threads, float delta) {
    if (debug) {
        printf("Starting solver with:\n"
               "n = %i\n"
               "iterations = %i\n"
               "threads = %i\n"
               "delta = %f\n",
               n, iterations, threads, delta);
    }

    // Allocate some buffers to calculate the solution in
    double *curr = (double*)calloc(n * n * n, sizeof(double));
    double *next = (double*)calloc(n * n * n, sizeof(double));
    double *temp;

    // Ensure we haven't run out of memory
    if (curr == NULL || next == NULL) {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    WorkerArgs args;
    args.xstart = 0;
    args.xend = n - 1;
    args.ystart = 0;
    args.yend = n - 1;
    args.zstart = 0;
    args.zend = n - 1;
    args.n = n;
    args.source = source;
    args.delta = delta;

    // TODO: solve Poisson's equation for the given inputs
    for (int i = 0; i < iterations; i++)
    {
        args.curr = curr;
        args.next = next;

        equation(&args);
        temp = curr;
        curr = next;
        next = temp;
    }

    // Free one of the buffers and return the correct answer in the other.
    // The caller is now responsible for freeing the returned pointer.
    free(next);

    if (debug) {
        printf("Finished solving.\n");
    }

    return curr;
}

void* worker(void* pargs, int n, double *source, int iterations, int threads, float delta) {
    WorkerArgs* args = (WorkerArgs*)pargs; //dont think this one does anything, but ill leave it here

    poisson_neumann(n, source, iterations, threads, delta);
    return NULL;
}



int main(int argc, char **argv) {

    // default settings for solver
    int iterations = 10;
    int n = 5;
    int threads = 1;
    float delta = 1;

    // parse the command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
            printf("usage: poisson [-n size] [-i iterations] [-t threads] [--debug]\n");
            return EXIT_SUCCESS;
        }

        if (strcmp(argv[i], "-n") == 0) {
            if (i == argc - 1) {
                fprintf(stderr, "Error: expected size after -n!\n");
                return EXIT_FAILURE;
            }

            n = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "-i") == 0) {
            if (i == argc - 1) {
                fprintf(stderr, "Error: expected iterations after -i!\n");
                return EXIT_FAILURE;
            }

            iterations = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "-t") == 0) {
            if (i == argc - 1) {
                fprintf(stderr, "Error: expected threads after -t!\n");
                return EXIT_FAILURE;
            }

            threads = atoi(argv[++i]);
        }

        if (strcmp(argv[i], "--debug") == 0) {
            debug = true;
        }
    }

    // ensure we have an odd sized cube
    if (n % 2 == 0) {
        fprintf(stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    double *source = (double*)calloc(n * n * n, sizeof(double));
    if (source == NULL) {
        fprintf(stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    source[(n * n * n) / 2] = 1;

    // Calculate the resulting field with Neumann conditions
    double *result = poisson_neumann(n, source, iterations, threads, delta);

    // Print out the middle slice of the cube for validation
        for (int x = 0; x < n; ++x) {
            for (int y = 0; y < n; ++y) {
                printf("%0.5f ", result[((n/2) * n + y) * n + x]);
            }
            printf("\n");
        }

    free(source);
    free(result);

    return EXIT_SUCCESS;
}
