#include "build.cpp"

#include <TApplication.h>
#include <sstream>
#include <iostream>

int main( int argc, const char* argv[] )
{
    TApplication * myapp = new TApplication("myapp", 0, 0);
    
    build();
    
    return EXIT_SUCCESS;
}
