import json


def trim(polynomial):
    while len(polynomial) > 1 and polynomial[-1] == 0:
        polynomial.pop()
    return polynomial


def add(left, right):
    return trim([(left[index] if index < len(left) else 0)
                 + (right[index] if index < len(right) else 0)
                 for index in range(max(len(left), len(right)))])


def scale(polynomial, factor):
    return trim([factor*value for value in polynomial])


def multiply(left, right):
    result = [0]*(len(left)+len(right)-1)
    for first, coefficient in enumerate(left):
        for second, value in enumerate(right):
            result[first+second] += coefficient*value
    return trim(result)


def power(polynomial, exponent):
    result = [1]
    for degree in range(exponent):
        result = multiply(result, polynomial)
    return result


def derivative(polynomial):
    return trim([degree*polynomial[degree] for degree in range(1,len(polynomial))] or [0])


def compose(polynomial, argument):
    result = [0]
    for coefficient in reversed(polynomial):
        result = add(multiply(result,argument),[coefficient])
    return result


RECURRENCES = {
    'A386830': {'r':2, 'q':3, 'numerator':[1], 'correction':[0,-30],
        'coefficients':[[0,-24,16],[-480,1242,-648],[3240,-7290,3645]]},
    'A386831': {'r':3, 'q':3, 'numerator':[1],
        'correction':[0,338670,-181329840,604049400],
        'coefficients':[[0,37488,-168696,168696],
            [-822702720,1658851582,-1040550807,203878391],
            scale([777640080,-1931547954,1277608165,-252935041],27),
            scale([523823440,-614869982,259381745,-38797609],3645),
            scale(multiply(multiply([-7,2],[-13,4]),[-11,4]),-16309333800)]},
    'A386766': {'r':1, 'q':2, 'numerator':[1,-5],
        'correction':[0,-75,300],
        'coefficients':[[0,9],[-30,-108],[2960,-1360],[4000,-1600]]}
}


def polynomial_identity(record):
    slope, weight = record['r'], record['q']
    first, second, base = [1,-weight], [1,-5*(slope+1)], [1,-5]
    coordinate = multiply([0,1],power(base,slope))
    order = max(len(polynomial)-1 for polynomial in record['coefficients'])
    numerators = [record['numerator']]
    for degree in range(order):
        inner = add(multiply(multiply(first,second),derivative(numerators[-1])),
                    multiply(add(scale(second,(degree+1)*weight),
                                 scale(first,(2*degree+1)*5*(slope+1))),numerators[-1]))
        numerators.append(multiply(multiply([0,1],base),inner))
    residual = [0]
    for shift, polynomial in enumerate(record['coefficients']):
        translated = compose(polynomial,[shift,1])
        for degree, coefficient in enumerate(translated):
            term = multiply(power(coordinate,shift),numerators[degree])
            term = multiply(term,power(first,order-degree))
            term = multiply(term,power(second,2*(order-degree)))
            residual = add(residual,scale(term,coefficient))
    correction = multiply(compose(record['correction'],coordinate),
                          multiply(power(first,order+1),power(second,2*order+1)))
    assert add(residual,scale(correction,-1)) == [0]
    return {'polynomial_numerator':[0], 'theta_order':order,
            'recurrence_onset':len(record['coefficients'])-1,
            'correction':record['correction']}


def reversion_identity():
    first, base, second = [1,-3],[1,-5],[1,-10]
    coordinate = multiply([0,1],base)
    theta_numerator = multiply([0,1],add(scale(base,3),scale(first,5)))
    term = multiply(compose([2,-31,-180],coordinate),theta_numerator)
    term = add(term,multiply(compose([2,-2,-270],coordinate),multiply(second,first)))
    correction = multiply(compose([2,30],coordinate),
                          multiply(second,multiply(power(first,2),base)))
    assert add(term,scale(correction,-1)) == [0]
    return {'polynomial_numerator':[0], 'recurrence_onset':2}


def main():
    identities = {name:polynomial_identity(record) for name,record in RECURRENCES.items()}
    identities['A386769'] = reversion_identity()
    print(json.dumps({'status':'pass', 'arithmetic':'integer polynomial identities',
                      'identities':identities},indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
