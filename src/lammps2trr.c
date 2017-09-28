#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <stdbool.h>
#include <string.h>
#include <argp.h>
#include <cblas.h>
#include "gromacs/fileio/gmxfio.h"
#include "gromacs/fileio/trrio.h"
#include "verbPrintf.c"

// CLI stuff
const char* argp_program_version = "lammps2trr";
const char* argp_program_bug_address = "<bernhardt@cpc.tu-darmstadt.de>";
static char doc[] = "\
lammps2trr -- convert lammps dump with velocities to trr file\n\
\n\
Does assume orthorombic box!\n\
Puts step as time in trr file!\n\
Does assume 'real' lammps units e.g. Angstrom and Angstrom/fs and converts them to gromacs units e.g. nm and nm/ps!\n";
static char args_doc[] = "";

static struct argp_option options[] = {
    {"verbose",  'v', 0,      0,  "Produce verbose output" },
    {"lammpstrj", 'f', "FILE", 0,  "Input lammps trajectory file with columns xu yu zu vx vy vz (default: traj.dump)"},
    {"trr", 'o', "FILE", 0,  "Output trr trajectory file (default: traj.trr)"},
    {"dt", 'd', "float", 0,  "Timestep in the lammps file, for the trr file (default: 0.001 ps)"},
    { 0 }
};

struct arguments
{
    bool verbosity;
    char *infile;
    char *outfile;
    float dt;
};

static error_t parse_opt (int key, char *arg, struct argp_state *state)
{
    /* Get the input argument from argp_parse, which we
       know is a pointer to our arguments structure. */
    struct arguments *arguments = state->input;

    switch (key)
    {
        case 'v':
            arguments->verbosity = true;
            break;
        case 'f':
            arguments->infile = arg;
            break;
        case 'o':
            arguments->outfile = arg;
            break;
        case 'd':
            arguments->dt = atof(arg);
            break;

        case ARGP_KEY_ARG:
            if (state->arg_num >= 0)
                /* Too many arguments. */
                argp_usage (state);

            break;

        case ARGP_KEY_END:
            if (state->arg_num < 0)
                /* Not enough arguments. */
                argp_usage (state);
            break;

        default:
            return ARGP_ERR_UNKNOWN;
    }
    return 0;
}

/* Our argp parser. */
static struct argp argp = {options, parse_opt, args_doc, doc};

int main( int argc, char *argv[] )
{
    // command line arguments
    struct arguments arguments;

    // Default values.
    arguments.verbosity = false;
    arguments.infile = "traj.dump";
    arguments.outfile = "traj.trr";
    arguments.dt = 0.001;

    // parse command line arguments
    argp_parse(&argp, argc, argv, 0, 0, &arguments);
    bool verbosity = arguments.verbosity;

    // variables
    int index_xuyuzuvxvyvz[6];
    const char* strings_xuyuzuvxvyvz[] = {"xu", "yu", "zu", "vx", "vy", "vz"};
    float lambda = 0.0;
    int natoms = 1;
    int step = 0;
    rvec* x = malloc(3 * natoms * sizeof(real));
    rvec* v = malloc(3 * natoms * sizeof(real));
    float lammps_box[6] = {0, 0, 0, 0, 0, 0}; // xmin, xmax, ymin, ...
    rvec box[3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

    // open lammps input file
    verbPrintf(verbosity, "opening file %s\n", arguments.infile);
    FILE* fp = fopen(arguments.infile, "r");
    if (fp == NULL)
    {
        perror("Opening lammps file for reading: ");
        return 1;
    }

    // parse first frame of lammps dump for natoms and position (column index) of xu..vz
    verbPrintf(verbosity, "read first frame from lammps file\n");
    char line[256];
    bool got_natoms = false;
    bool got_coordpos = false;
    while (fgets(line, sizeof(line), fp))
    {
        //printf("%s", line);

        if (strncmp(line, "ITEM", 4) == 0)
        {
            if (strncmp(line, "ITEM: NUMBER OF ATOMS", 21) == 0)
            {
                fgets(line, sizeof(line), fp);
                sscanf(line, "%d", &natoms);
                x = malloc(3 * natoms * sizeof(real));
                v = malloc(3 * natoms * sizeof(real));
                got_natoms = true;
            }
            else if (strncmp(line, "ITEM: ATOMS", 11) == 0)
            {
                //loop over tokens of line
                verbPrintf(verbosity, "column indices found:\n");
                char* token;
                int column = 0;
                token = strtok(line, " ");
                while (token != NULL)
                {
                    // loop over stringst to compare with
                    for (int i=0; i<6; i++)
                    {
                        // compare token with string
                        if (strncmp(token, strings_xuyuzuvxvyvz[i], 2) == 0)
                        {
                            // set index to column
                            index_xuyuzuvxvyvz[i] = column - 2;
                            verbPrintf(verbosity, "%s %d\n", strings_xuyuzuvxvyvz[i], index_xuyuzuvxvyvz[i]);
                        }
                    }

                    column++;
                    token = strtok(NULL, " ");
                }
                got_coordpos = true;
            }
        }

        if (got_natoms && got_coordpos) break;
    }

    // go back to file start
    rewind(fp);

    // open trr output file
    verbPrintf(verbosity, "starting writing file %s\n", arguments.outfile);
    t_fileio* trj_out = gmx_trr_open(arguments.outfile, "w"); /* should check the result */

    // parse lammps dump
    verbPrintf(verbosity, "read frames from lammps file\n");
    while (fgets(line, sizeof(line), fp))
    {
        if (strncmp(line, "ITEM", 4) == 0)
        {
            if (strncmp(line, "ITEM: NUMBER OF ATOMS", 21) == 0)
            {
                // natoms already parsed in first frame
                fgets(line, sizeof(line), fp);
            }
            else if (strncmp(line, "ITEM: TIMESTEP", 14) == 0)
            {
                fgets(line, sizeof(line), fp);
                sscanf(line, "%d", &step);
            }
            else if (strncmp(line, "ITEM: BOX BOUNDS pp pp pp", 25) == 0)
            {
                fgets(line, sizeof(line), fp);
                sscanf(line, "%f %f", &lammps_box[0], &lammps_box[1]);

                fgets(line, sizeof(line), fp);
                sscanf(line, "%f %f", &lammps_box[2], &lammps_box[3]);

                fgets(line, sizeof(line), fp);
                sscanf(line, "%f %f", &lammps_box[4], &lammps_box[5]);

                // convert box (orthorombic only) + unit conversion
                box[0][0] = (lammps_box[1] - lammps_box[0]) / 10;
                box[1][1] = (lammps_box[3] - lammps_box[2]) / 10;
                box[2][2] = (lammps_box[5] - lammps_box[4]) / 10;
            }
            else if (strncmp(line, "ITEM: ATOMS", 11) == 0)
            {
                for (int j=0; j<natoms; j++)
                {
                    fgets(line, sizeof(line), fp);
                    char* token;
                    int column = 0;
                    token = strtok(line, " ");
                    while (token != NULL)
                    {
                        // loop over strings to compare with
                        for (int i=0; i<6; i++)
                        {
                            // if column is interesting
                            if (column == index_xuyuzuvxvyvz[i])
                            {
                                // token to array + box shift to (0 0 0) + unit conversion
                                switch(i) {
                                    case 0: x[j][0] = (atof(token) - lammps_box[0]) / 10; break;
                                    case 1: x[j][1] = (atof(token) - lammps_box[2]) / 10; break;
                                    case 2: x[j][2] = (atof(token) - lammps_box[4]) / 10; break;
                                    case 3: v[j][0] = atof(token) * 100; break;
                                    case 4: v[j][1] = atof(token) * 100; break;
                                    case 5: v[j][2] = atof(token) * 100; break;
                                }
                            }
                        }

                        column++;
                        token = strtok(NULL, " ");
                    }
                    got_coordpos = true;
                }

                // write to trr trajectory
                gmx_trr_write_frame(trj_out, step, step * arguments.dt, lambda, box, natoms, x, v, NULL);
            }
        }
        else
        {
            printf("THIS SHOULD NOT HAPPEN! Unexpected data in lammps dump?\n");
        }


    }
    fclose(fp);





    gmx_trr_close(trj_out);

    return 0;
}
