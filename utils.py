def calculate_greeks(pricing_function, S0, K, r, T, sigma, num_simulations, option_type, **kwargs):
    epsilon = 1e-4
    price = pricing_function(S0, K, r, T, sigma, num_simulations, option_type, **kwargs)

    price_up = pricing_function(S0 + epsilon, K, r, T, sigma, num_simulations, option_type, **kwargs)
    price_down = pricing_function(S0 - epsilon, K, r, T, sigma, num_simulations, option_type, **kwargs)
    delta = (price_up - price_down) / (2 * epsilon)
    gamma = (price_up - 2 * price + price_down) / (epsilon ** 2)

    price_shorter = pricing_function(S0, K, r, T - epsilon, sigma, num_simulations, option_type, **kwargs)
    theta = (price_shorter - price) / epsilon

    price_vol_up = pricing_function(S0, K, r, T, sigma + epsilon, num_simulations, option_type, **kwargs)
    vega = (price_vol_up - price) / epsilon

    price_rate_up = pricing_function(S0, K, r + epsilon, T, sigma, num_simulations, option_type, **kwargs)
    rho = (price_rate_up - price) / epsilon

    return delta, gamma, theta, vega, rho
