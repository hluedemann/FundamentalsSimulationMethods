/** Fundamentals of Simulation Methods
*
*   Author: H.Lüdemann, F.Walberg, D.Wolf
*   Exercise: Sheet 1 part 2
*   Topic: Near-cancellation of numbers
*   Due date: 23.10.2019
**/


#include <iostream>
#include <cmath>


/**
 * @brief Function to compute
 */
double f(double x)
{
    return (x + std::exp(-x) - 1) / (x * x);
}

/**
 * @brief Function to compute with adjustments for small values.
 * This verision of the function calculates the taylor expansion if the value
 * gets small.
 */
double adjustedF(double x)
{
    double minValue = 5.0*10e-6 ;

    if(x > minValue)
    {
        return (x + std::exp(-x) - 1) / (x * x);
    }
    else
    {
        return 0.5 - x / 6.0 + x * x / 24.0;
    }
}

/**
 * @brief Promt the user to input a value.
 * @return The value given by the user.
 */
double getUserInput()
{
    std::string input;
    double x;
    std::cout << "Enter a value for x." << std::endl;
    std::cin >> input;

    try
    {
        x = std::stod(input);
    }
    catch(...)
    {
        std::cerr << "Invalid value entered please enter again." << std::endl;;
        x = getUserInput();
    }
    
    return x;
}

/**
 * @brief Evaluate the function for the given user input and print the result.
 */
void calculateUserInput()
{
    double x = getUserInput();

    // Change the function to f(x) if you want to see the wrong results for small values.
    double result = adjustedF(x);

    std::cout << "Result:" << std::endl;
    std::cout << "f(" << x << ") = " << result << std::endl;
}


/**
 * @brief Evaluate the function f(x) for small values.
 * This function is used to determin the lower limit for which the function f(x)
 * yields correct results.
 */
void printResultForSmallValues()
{
    std::cout << "#x\ty" << std::endl;

    double minValue = 10e-10;

    for(double i = 1.0; i > minValue; i /= 1.5)
    {
        std::cout << i << "\t" << adjustedF(i) << std::endl;
    }
}


int main(int argc, char** argv)
{
    int id = 0;
    if(argc > 1)
    {
        id = atoi(argv[1]);
    }
    // Adjust this function to evaluate f(x) or adjustedF(x)
    if(id==0)
    {
        calculateUserInput();
    }
    // Adjust this function to ese f(x) or adjustedF(x).
    else if(id == 1)
    {
        printResultForSmallValues();
    }
    else
    {
        printf("Id 0: get user input\nId 1: small value breakdown\n");
    }
    
}

