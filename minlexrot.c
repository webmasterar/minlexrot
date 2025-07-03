/*
* Copyright (c) 2025 Ahmad Retha; MIT License.
*/

#include <argp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "generalised_min_lex_rot.h"

static char doc[] =
    "minlexrot -- A generalised implementation of the Lexicographically Minimal String Rotation using algorithm BW or KR";

static char args_doc[] =
    "alphabet\t-a\tThe alphabet used for the string, e.g. 'ACGT'.\n"
    "algorithm\t-g\tThe algorithm to use for calculation (BW=Bitwise or KR=Karp-Rabin).\n"
    "direction\t-d\tThe direction of rotation (FW=Forward or BW=Backward). Defaults to Backward.\n"
    "text\t-t\tThe input string for which to calculate the lexicographically minimal rotation.";

static struct argp_option options[] = {
    {"alphabet", 'a', "ALPHABET", 0, "The alphabet used for the string, relevant for bitwise operations."},
    {"algorithm", 'g', "ALGORITHM", 0, "The algorithm to use for calculation (BW=Bitwise or KR=Karp-Rabin)."},
    {"direction", 'd', "DIRECTION", 0, "The direction of rotation (FW=Forward or BW=Backward). Defaults to Backward."},
    {"text", 't', "TEXT", 0, "The input string for which to calculate the lexicographically minimal rotation."},
    { 0 }
};

struct arguments
{
    char *alphabet;
    char *algorithm;
    char *direction;
    char *text;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'a':
            arguments->alphabet = arg;
            break;
        case 'g':
            if (strcmp(arg, "BW") != 0 && strcmp(arg, "KR") != 0) {
                argp_error(state, "Algorithm must be 'BW' or 'KR'");
            }
            arguments->algorithm = arg;
            break;
        case 'd':
            if (strcmp(arg, "FW") != 0 && strcmp(arg, "BW") != 0) {
                argp_error(state, "Direction must be 'FW' or 'BW'");
            }
            arguments->direction = arg;
            break;
        case 't':
            arguments->text = arg;
            break;

        case ARGP_KEY_ARG:
            return ARGP_ERR_UNKNOWN; // No non-option arguments expected
        case ARGP_KEY_END:
            // Check for required arguments
            if (!arguments->alphabet)
            argp_error (state, "Required argument --alphabet is missing");
            if (!arguments->algorithm)
            argp_error (state, "Required argument --algorithm is missing");
            if (!arguments->direction)
            argp_error (state, "Required argument --direction is missing");
            if (!arguments->text)
            argp_error (state, "Required argument --text is missing");
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }

    return 0;
}

static struct argp argp = { options, parse_opt, args_doc, doc };

int main (int argc, char **argv)
{
    struct arguments arguments;

    arguments.alphabet = NULL;
    arguments.algorithm = NULL;
    arguments.direction = NULL;
    arguments.text = NULL;

    argp_parse(&argp, argc, argv, 0, 0, &arguments);

    uint32_t n = strlen(arguments.text);
    uint32_t a = strlen(arguments.alphabet);
    uint32_t algo = (strcmp(arguments.algorithm, "BW") == 0) ? ALGO_BW : ALGO_KR;
    uint32_t roll_dir = (strcmp(arguments.direction, "BW") == 0) ? ROLL_DIR_BACKWARD : ROLL_DIR_FORWARD;

    uint32_t rot = min_lex_rot(arguments.text, n, arguments.alphabet, a, algo, roll_dir);

    char* T = (char*) malloc(n + 1);
    uint32_t i;
    for (i = 0; i < n; i++) {
        T[i] = arguments.text[(rot + i) % n];
    }
    T[n] = '\n';

    printf("%d %s\n", rot, T);

    free(T);

    return 0;
}
