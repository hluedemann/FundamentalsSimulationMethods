/** Fundamentals of Simulation Methods
*
*   Author: Hauke Lüdemann
*   Exercise: Sheet 1 part 1
*   Topic: Pitfalls of integer and floating point arithmetic
*   Due date: 23.10.2019
**/


#include <iostream>


/**
 * @brief This is a simple program to demonstrate some pitfalls whenn calculating with ints and floats.
 */
int main()
{
    /// Number 1.1
    {
        int i = 7;
        float y = 2 * (i / 2);
        float z = 2 * (i / 2.);

        printf("Part 1.1:\n");
        printf("%e %e \n", y, z);
    }

    /// Number 1.2
    {
        double a = 1.0e17;
        double b = -1.0e17;
        double c = 1.0;
        double x = (a + b) + c;
        double y = a + (b + c);

        printf("Part 1.2:\n");
        printf("%e %e \n", x, y);
    }

    /// Number 1.3
    {
        float x = 1e20;
        float y;
        y = x * x;

        printf("Part 1.3:\n");
        printf("%e %e\n", x, x / y);
    }
    return 0;
}
