#include "main.h"

int main(void) {
    // Initialize system
    System system;
    init_system(&system, SYSTEM_METHOD, SYSTEM_YEAR, SYSTEM_SIZE, SYSTEM_STEP, SYSTEM_TIME);

    // Run the simulation
    run(&system);

    // Create file for results
    FILE *json = fopen(NAME, "wt");
    if (json == NULL) {
        printf("ERROR: Could not open file!\n");
        return EXIT_FAILURE;
    }

    // Write data to file
    write_data(&system, json);

    return EXIT_SUCCESS;
}
