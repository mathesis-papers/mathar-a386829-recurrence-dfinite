import json
from math import comb
from fractions import Fraction
from verify_polynomials import RECURRENCES


def multiply(left, right, length):
    result = [0] * length
    for degree in range(length):
        result[degree] = sum(left[index] * right[degree - index]
                             for index in range(degree + 1)
                             if index < len(left) and degree - index < len(right))
    return result


def inverse_unit(series, length):
    assert series[0] == 1
    result = [1] + [0] * (length - 1)
    for degree in range(1, length):
        result[degree] = -sum(series[index] * result[degree - index]
                              for index in range(1, min(degree + 1, len(series))))
    return result


def defining_sum(index, slope=1, first=3, second=2):
    return sum(first**part * second**(index-part)
               * comb((slope+1)*index+1, part)
               * comb((slope+1)*index-part, index-part)
               for part in range(index+1))


def auxiliary_sum(index, slope=1, total=5, first=3):
    return sum(total**part * first**(index-part) * comb(slope*index+part, part)
               for part in range(index+1))


def alternative_sum(index, slope=1, total=5, second=2):
    return sum(total**part * (-second)**(index-part)
               * comb((slope+1)*index+1, part) for part in range(index+1))


def implicit_series(slope, length):
    result = [0] * length
    for degree in range(1, length):
        numerator = 5**(degree-1) * comb((slope+1)*degree-2, degree-1)
        assert numerator % degree == 0
        result[degree] = numerator // degree
    return result


def family_series(slope, length):
    root = implicit_series(slope, length)
    first = [1] + [-3*value for value in root[1:]]
    second = [1] + [-5*(slope+1)*value for value in root[1:]]
    return inverse_unit(multiply(first, second, length), length)


def shifted_sum(index, slope=1, total=5, first=3):
    if index == 0:
        return 1
    return sum(total**part*first**(index-part)*comb(slope*index+part-1,part)
               for part in range(index+1))


def reversion_sum(index, slope=1, first=3, second=2):
    numerator = sum(first**part*second**(index-part)
                    *comb((slope+1)*(index+1),part)
                    *comb((slope+1)*index+slope-part-1,index-part)
                    for part in range(index+1))
    assert numerator % (index+1) == 0
    return numerator//(index+1)


def expanded_checks(length):
    ranges = {}
    for slope in (1,2,3):
        for first,second in ((3,2),(2,3),(1,2)):
            total = first+second
            direct = [defining_sum(index,slope,first,second) for index in range(length)]
            shifted = [shifted_sum(index,slope,total,first) for index in range(length)]
            reverted = [reversion_sum(index,slope,first,second) for index in range(length)]
            root = [0]+[total**(degree-1)*comb((slope+1)*degree-2,degree-1)//degree
                        for degree in range(1,length)]
            first_factor = [1]+[-first*value for value in root[1:]]
            second_factor = [1]+[-(slope+1)*total*value for value in root[1:]]
            assert direct == inverse_unit(multiply(first_factor,second_factor,length),length)
            numerator = [1]+[-total*value for value in root[1:]]
            assert shifted == multiply(numerator,direct,length)
            for index in range(length):
                assert direct[index] == auxiliary_sum(index,slope,total,first)
                assert direct[index] == alternative_sum(index,slope,total,second)
                assert shifted[index] == sum(total**part*(-second)**(index-part)
                    *comb((slope+1)*index,part) for part in range(index+1))
            for index in range(1,length):
                assert index*reverted[index] == sum(shifted[part]*reverted[index-part]
                                                       for part in range(1,index+1))
                ratio = Fraction(total,first)
                hypergeom = (1-ratio)**(-slope*index)*sum(
                    Fraction(index,index+part)*(-1)**part*comb(slope*index-1,part)*ratio**part
                    for part in range(slope*index))
                expression = first**index*(1-ratio)**(-slope*index)
                expression -= total**index*comb((slope+1)*index-1,index)*(hypergeom-1)
                assert expression == shifted[index]
            ranges[f'r={slope},u={first},v={second}'] = [0,length-1]
    for name,record in RECURRENCES.items():
        values = [defining_sum(index,record['r']) if name != 'A386766'
                  else shifted_sum(index,first=2) for index in range(length)]
        onset = len(record['coefficients'])-1
        for index in range(onset,length):
            residual = sum(sum(coefficient*index**degree
                               for degree,coefficient in enumerate(polynomial))*values[index-shift]
                           for shift,polynomial in enumerate(record['coefficients']))
            assert residual == 0
    reverted = [reversion_sum(index) for index in range(length)]
    for index in range(2,length):
        assert 2*(index+1)*reverted[index]+(-31*index+29)*reverted[index-1]+90*(-2*index+1)*reverted[index-2] == 0
    return ranges


def main():
    length = 81
    sequence = [defining_sum(index) for index in range(length)]
    assert sequence[:6] == [1, 13, 204, 3457, 61006, 1103598]
    central = [5**index * comb(2*index, index) for index in range(length)]
    for index in range(length):
        assert sequence[index] == auxiliary_sum(index) == alternative_sum(index)
        assert 4*sequence[index] + (18*sequence[index-1] if index else 0) == 7*central[index] - (3 if index == 0 else 0)
    for index in range(2, length):
        assert 2*index*sequence[index] + (20-31*index)*sequence[index-1] + 90*(1-2*index)*sequence[index-2] == 0
    for degree in range(length-1):
        residual = 2*(degree+1)*sequence[degree+1] - (31*degree+11)*sequence[degree]
        if degree:
            residual -= (180*(degree-1)+270)*sequence[degree-1]
        assert residual == (15 if degree == 0 else 0)
    for degree in range(length-2):
        residual = 2*(degree+2)*(degree+1)*sequence[degree+2] - (degree+1)*(31*degree+42)*sequence[degree+1]
        residual -= (180*degree*(degree-1)+630*degree+270)*sequence[degree]
        assert residual == 0
    numerator = multiply([4, 18], sequence, length)
    numerator[0] += 3
    algebraic = multiply([1, -20], multiply(numerator, numerator, length), length)
    assert algebraic == [49] + [0]*(length-1)
    sibling_length = 31
    swapped = [defining_sum(index, first=2, second=3) for index in range(sibling_length)]
    for index in range(sibling_length):
        assert swapped[index] == auxiliary_sum(index, first=2)
        assert swapped[index] == alternative_sum(index, second=3)
        assert 3*swapped[index] + (4*swapped[index-1] if index else 0) == 4*central[index] - (1 if index == 0 else 0)
    for index in range(2, sibling_length):
        assert 3*index*swapped[index] + (30-56*index)*swapped[index-1] - 40*(2*index-1)*swapped[index-2] == 0
    shifted = [1] + [sum(5**part * 3**(index-part) * comb(index+part-1, part)
                         for part in range(index+1)) for index in range(1, sibling_length)]
    assert shifted[:5] == [1, 8, 114, 1862, 32246]
    for index in range(sibling_length):
        assert 2*shifted[index] + (9*shifted[index-1] if index else 0) == central[index] + (15*central[index-1] if index else 1)
    for slope in [1, 2, 3]:
        root = implicit_series(slope, sibling_length)
        factor = [1] + [-5*value for value in root[1:]]
        implicit = root
        for exponent in range(slope):
            implicit = multiply(implicit, factor, sibling_length)
        assert implicit == [0, 1] + [0]*(sibling_length-2)
        formula = family_series(slope, sibling_length)
        assert formula == [defining_sum(index, slope) for index in range(sibling_length)]
        assert formula == [auxiliary_sum(index, slope) for index in range(sibling_length)]
        assert formula == [alternative_sum(index, slope) for index in range(sibling_length)]
    print(json.dumps({
        'status': 'pass',
        'arithmetic': 'exact integers and rational numbers',
        'parameter_families': expanded_checks(31),
        'added_recurrence_ranges': {'A386830':[2,30], 'A386831':[4,30],
                                    'A386766':[3,30], 'A386769':[2,30]},
        'exponential_reversion_and_hypergeometric_ranges':[1,30],
        'A386829': {'coefficient_range': [0, 80], 'recurrence_range': [2, 80],
                    'inhomogeneous_equation_degrees': [0, 79],
                    'homogeneous_equation_degrees': [0, 78],
                    'algebraic_equation_degrees': [0, 80], 'initial_terms': sequence[:6]},
        'A385669': {'coefficient_range': [0, 30], 'recurrence_range': [2, 30]},
        'A386763': {'coefficient_range': [0, 30], 'initial_terms': shifted[:5]},
        'A386830': {'coefficient_range': [0, 30], 'implicit_root_equation': 'T*(1-5*T)^2=z'},
        'A386831': {'coefficient_range': [0, 30], 'implicit_root_equation': 'T*(1-5*T)^3=z'},
        'scope': 'Finite coefficient identities; proofs for all indices are in note.pdf.'
    }, indent=2))


if __name__ == '__main__':
    main()
