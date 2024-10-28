class Model:
    def __init__(self):
        raise NotImplementedError

    def get_zcb_mse(self):
        raise NotImplementedError

    def get_caplet_mse(self):
        raise NotImplementedError

    def calibrate_model(self):
        raise NotImplementedError

    def price_calplet(self):
        raise NotImplementedError

    def price_cap(self):
        raise NotImplementedError

    def price_swaption(self):
        raise NotImplementedError


class Curve:
    # For interest rate and volatility curves
    def __init__(self):
        raise NotImplementedError

    def interp(self):
        raise NotImplementedError


class Price_curve:
    def __init__(self):
        raise NotImplementedError
