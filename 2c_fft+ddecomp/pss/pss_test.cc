#include "pss_read.hh"

int main() {

    // Set up the perfect squared square class
    pss_read p(112,"pss.txt");

    // Print an ASCII art represention of the gridpoint types 
    p.ascii_art();
    
    // Analyze the gridpoints in the glueing region between squares
    p.analyze_glue();
}
