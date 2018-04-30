#include "car_model.hh"

int main() {

    // Initialize the car model with 60 cars. The length of the domain
    // will be set to 90.
    car_model cm(60);

    // Simulate to t=100, and output 200 snapshots
    cm.solve_fixed(100.,200,true);
}
