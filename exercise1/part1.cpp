#include <iostream>

int main()
{
    /// Number 1.1
    {
        int i = 7;
        float y = 2 * (i / 2);      // -> Integer devicion cuts of the floating part
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
        y = x * x;      // -> Overflow, hence devision by zero

        printf("Part 1.3:\n");
        printf("%e %e\n", x, x / y);
    }
    return 0;
}