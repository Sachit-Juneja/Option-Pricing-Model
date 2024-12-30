import numpy as np

def european_option_pricing(S0, K, r, T, sigma, num_simulations, option_type):
    Z = np.random.standard_normal(num_simulations)
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    if option_type == "call":
        payoffs = np.maximum(ST - K, 0)
    else:
        payoffs = np.maximum(K - ST, 0)
    return np.mean(payoffs) * np.exp(-r * T)

def asian_option_pricing(S0, K, r, T, sigma, num_simulations, option_type):
    dt = T / 100
    ST_path = S0 * np.cumprod(np.exp((r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * np.random.randn(100, num_simulations)), axis=0)
    avg_price = np.mean(ST_path, axis=0)
    if option_type == "call":
        payoffs = np.maximum(avg_price - K, 0)
    else:
        payoffs = np.maximum(K - avg_price, 0)
    return np.mean(payoffs) * np.exp(-r * T)

def barrier_option_pricing(S0, K, r, T, sigma, num_simulations, option_type, barrier, barrier_type):
    Z = np.random.standard_normal(num_simulations)
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    if barrier_type == "up-and-out":
        payoffs = np.maximum(ST - K, 0) * (ST < barrier)
    elif barrier_type == "down-and-out":
        payoffs = np.maximum(ST - K, 0) * (ST > barrier)
    elif barrier_type == "up-and-in":
        payoffs = np.maximum(ST - K, 0) * (ST >= barrier)
    elif barrier_type == "down-and-in":
        payoffs = np.maximum(ST - K, 0) * (ST <= barrier)
    else:
        raise ValueError("Invalid barrier type.")
    return np.mean(payoffs) * np.exp(-r * T)
