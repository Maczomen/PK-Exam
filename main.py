import numpy as np
from numpy.polynomial import Polynomial as P


def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    else:
        (d, x, y) = extended_gcd(b, a % b)
        return d, y, x - (a // b) * y


def modular_multiplication_inverse(num, mod):
    gcd, x, _ = extended_gcd(num, mod)
    if gcd == 1:
        return x
    else:
        raise ValueError("The modular inverse does not exist.")


def modular_addition_inverse(num, mod):
    return (mod - num) % mod


def hex_to_table(num):
    result = []
    while num != 0:
        result.append(num % 2)
        num //= 2
    return result


def extended_gcd_polynomial(a, b, order):
    if b == P([0]):
        return a, P([1]), P([0])
    else:
        # a_new = polynomial_mod_int((a % b), order)
        q, a_new = polynomial_floor_division_mod(a, b, order)
        d, s, t = extended_gcd_polynomial(b, a_new, order)
        # q3 = a % b
        # q = a // b
        # q2, q4 = polynomial_floor_division_mod(a, b, order)
        return d, t, s - q * t


def polynomial_floor_division_mod(pol1, pol2, mod):
    pol1_arr = list(reversed(pol1.coef))
    pol2_arr = list(reversed(pol2.coef))

    result = [0 for _ in range(len(pol1_arr) - len(pol2_arr) + 1)]
    for i in range(len(result)):
        if pol1_arr[0] == 0:
            pol1_arr.pop(0)
            continue

        divider = (modular_multiplication_inverse(pol2_arr[0], mod) * pol1_arr[0]) % mod
        for j, element in enumerate(pol2_arr):
            pol1_arr[j] -= divider * pol2_arr[j]
            pol1_arr[j] %= mod

        result[i] = divider
        if pol1_arr[0] == 0:
            pol1_arr.pop(0)
            continue
        else:
            raise ValueError("Division error.")

    if len(pol1_arr) == 0:
        pol1_arr = [0]
    if len(result) == 0:
        result = [0]
    return P(list(reversed(result))).trim(), P(list(reversed(pol1_arr))).trim()


def find_inverse_polynomial(element, modulus, order):
    element = P(element).trim()
    modulus = P(modulus).trim()
    d, s, t = extended_gcd_polynomial(element, modulus, order)
    test = s
    if d != P([1]):
        test, _ = polynomial_floor_division_mod(s, d, order)
    return polynomial_mod_int(test, order)


def polynomial_mod_int(polynomial, num):
    return P(np.mod(polynomial.coef, num)).trim()


def str_shorter(polynomial):
    result_string = ""
    for i, coef in enumerate(polynomial.coef):
        if coef != 0:
            result_string += f"{coef} * x^{i}"
            if i != polynomial.degree():
                result_string += " + "

    return result_string


if __name__ == '__main__':
    num = 9
    mod = 889
    try:
        inverse = modular_multiplication_inverse(num, mod) % mod
        print("Element odwrotny dla liczby", num, "w grupie (Z_" + str(mod) + ", ∗_" + str(mod) + ") wynosi:", inverse)
        print("Test:", num, "*", inverse, "%", mod, "=", num * inverse % mod)  # should return 1
    except ValueError as e:
        print("Błąd:", str(e))

    element_from_hex = hex_to_table(0x13)  # dla Z_2 [X] element = [1, 1, 0, 0, 1] można zapisać jako 0x1E

    element = [1, 1, 0, 0, 1, 0, 1]  # x^0, x^1, x^2, ....
    modulus = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]  # x^0, x^1, x^2, ....
    # element = hex_to_table(0x1)  # x^0, x^1, x^2, ....
    # element = [2, 2]  # x^0, x^1, x^2, ....
    # modulus = hex_to_table(0x5)  # x^0, x^1, x^2, ....
    order = 2

    inverse = find_inverse_polynomial(element, modulus, order)
    print("Element odwrotny dla polynomialu ")
    print(str_shorter(P(element)))
    print("W ciele Z_" + str(order) + " [X]/<" + str_shorter(P(modulus)) + "> wynosi: ")
    print(str_shorter(inverse))

    print("Test: " + str_shorter(polynomial_mod_int(inverse * element % modulus, order)))  # should return 1.0 * x^0
