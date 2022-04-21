#include <stdio.h>
#include <stdlib.h>

typedef struct{
    int n[6];       // 24 bytes
    double a[6];    // 48 bytes
    int c;          //  4 bytes
    unsigned int x; //  4 bytes
    unsigned int y; //  4 bytes
} header;
// either 80 or 84


typedef struct{
    int npart[6];                    // number of particles of each type in this file
    double mass[6];                  /* mass of particles of each type. If 0, then the
                                                                    masses are explicitly stored in the mass-block of
                                                                    the snapshot file, otherwise they are omitted */
    double time;                     // Time of snapshot file
    double redshift;                 // Redshift of snapshot file
    int flag_sfr;                    // Simulation includes star formation
    int flag_feedback;           // Feedback included (obsolete)
    unsigned int npartTotal[6];  // Number of particles of each type in this snapshot.
    int flag_cooling;            // Cooling included
    int num_files;               // Number of files in multi-file snapshot
    double BoxSize;              // Box-size of simulation
    double Omega0;               // matter density in units of critical density
    double OmegaLambda;          // cosmological constant parameter
    double HubbleParam;          // Hubble parameter in units of 100 km/sec/Mpc
    int flag_stellarage;         // File contains formation times of star particles
    int flag_metals;             // File contains metallicity for gas and star particles
    unsigned int npartTotalHighWord[6];  /* High word of the total number of particles
                                           of each type */
    int  flag_entropy_instead_u; /* flags that IC-file contains entropy instead
                                  of u */
    int flag_RII;               //TypeII SN rate per year
    int flag_RIa;               //TypeIa SN rate per year
    //char flag_Carbon;
    //char flag_Nitrogen;
    //char flag_Oxygen;
    //char flag_Florine;
    //char flag_Neon;
    //char flag_Sodium;
    //char flag_Magnesium;
    //char flag_Aluminum;
    //char flag_Silicon;
    //char flag_Phosphorus;
    //char flag_Sulfur;
    //char flag_Chlorine;
    //char flag_Argon;
    //char flag_Potassium;
    //char flag_Calcium;
    //char flag_Scandium;
    //char flag_Titanium;
    //char flag_Vanadium;
    //char flag_Chromium;
    //char flag_Manganese;
    //char flag_Iron;
    //char flag_Cobalt;
    //char flag_Nickel;
    //char flag_Copper;
    //char flag_Zinc;
    //char fill[27];               // fills to 256 Bytes
} gadget2_header;


int main(){
    printf("sizeof(header) = %li\n", sizeof(header));
    printf("sizeof(gadget2_header) = %li\n", sizeof(gadget2_header));

    return 0;
}
