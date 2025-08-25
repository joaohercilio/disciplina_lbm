#include "Simulation.hpp"

int main(int argc, char* argv[]) 
{    
    auto simulation = Simulation::create();

    simulation->run();  
}