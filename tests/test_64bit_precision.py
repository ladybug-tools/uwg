import pytest
from decimal import Decimal

def test_64bit_precision():
    """
    Tests for max number of precision for 64bit python float
    64 bits = 53 bits (significand) + 8 bits (characteristic)
    Therefore max fraction size for 53 bits in base 2 is (2^-53 = 1.11e-16)
    So for numbers > 0. we have to subtract the 16 significant digits (significand)
    from the integer part when testing the numbers
    """

    """
    # Unblock this for visual test
    print "0.12345678901234561111", "-- original input"
    print Decimal.from_float(0.12345678901234561111)
    print "  1234567890123456 -- will be accurate to 16th digit"
    print Decimal.from_float(11.12345678901234561111)
    print "   12345678901234 -- will be accurate to 14th digit"
    print Decimal.from_float(11111.12345678901234561111)
    print "      12345678901 -- will be accurate to 11th digit"
    print '----'
    """

    # 16 digits number
    digit_19 = 0.1234567890123456444
    trunc_at_11 = digit_19*10**11 - int(digit_19*10**11) # arbritrary truncation
    # trunc_at_11: 0.23456444
    #print Decimal.from_float(trunc_at_11)
    # != at 17+1 digits (give one additional digit for tolerance)
    assert trunc_at_11 != pytest.approx(0.23456444, abs=1e-7)
    # == at 16 digits
    assert trunc_at_11 == pytest.approx(0.23456444, abs=1e-5)

    # Equal at 11 digits
    digit_11 = 11111.12345678901444
    trunc_at_5 = digit_11*10**5 - int(digit_11*10**5) # arbritrary truncation
    # trunc_at_5: 0.678901444
    #print Decimal.from_float(trunc_at_5)
    # != at 17+1 digits (give one additional digit for tolerance)
    assert trunc_at_5 != pytest.approx(0.678901444, abs=1e-8)
    # == at 16 digits
    assert trunc_at_5 == pytest.approx(0.678901444, abs=1e-6)


if __name__ == "__main__":
    test_64bit_precision()
