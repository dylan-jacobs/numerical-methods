clc;
clear variables;

syms A B C deltaX 

eqn1 = A + B + C == 0;
eqn2 = (B * deltaX) + (C * 2 * deltaX) == 0;
eqn3 = (B * (deltaX ^ 2) / 2) + (2 * C * (deltaX ^ 2)) - 1 == 0;

solution = solve([eqn1, eqn2, eqn3], [A, B, C])
