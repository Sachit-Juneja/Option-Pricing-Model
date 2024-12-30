import numpy as np

def monte_carlo_option_pricing_with_full_greeks(S0, K, r, T, sigma, num_simulations, option_type="call"):
    """
    Monte Carlo simulation for European options pricing with full Greeks.

    Parameters:
    - S0: Initial stock price
    - K: Strike price
    - r: Risk-free interest rate
    - T: Time to maturity (in years)
    - sigma: Volatility of the underlying asset
    - num_simulations: Number of simulations
    - option_type: "call" or "put"

    Returns:
    - Estimated option price and Greeks (Delta, Gamma, Theta, Vega, Rho)
    """
    # Generate random standard normal numbers
    Z = np.random.standard_normal(num_simulations)
    
    # Simulate end prices for the stock
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    
    # Calculate payoffs for the option
    if option_type == "call":
        payoffs = np.maximum(ST - K, 0)
    elif option_type == "put":
        payoffs = np.maximum(K - ST, 0)
    else:
        raise ValueError("Invalid option type. Choose 'call' or 'put'.")
    
    # Discount payoffs back to present value
    discounted_payoffs = payoffs * np.exp(-r * T)
    
    # Estimate the option price
    option_price = np.mean(discounted_payoffs)
    
    # Delta (sensitivity to S0)
    epsilon = 1e-4 * S0
    S0_up = S0 + epsilon
    S0_down = S0 - epsilon
    ST_up = S0_up * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    ST_down = S0_down * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    
    if option_type == "call":
        payoff_up = np.maximum(ST_up - K, 0)
        payoff_down = np.maximum(ST_down - K, 0)
    elif option_type == "put":
        payoff_up = np.maximum(K - ST_up, 0)
        payoff_down = np.maximum(K - ST_down, 0)
    
    delta = (np.mean(payoff_up) - np.mean(payoff_down)) / (2 * epsilon) * np.exp(-r * T)
    
    # Gamma (sensitivity of Delta to changes in S0)
    gamma = (np.mean(payoff_up) - 2 * np.mean(payoffs) + np.mean(payoff_down)) / (epsilon ** 2) * np.exp(-r * T)
    
    # Theta (sensitivity to time to maturity)
    epsilon_t = 1e-4
    T_down = T - epsilon_t
    ST_t_down = S0 * np.exp((r - 0.5 * sigma**2) * T_down + sigma * np.sqrt(T_down) * Z)
    if option_type == "call":
        payoff_t_down = np.maximum(ST_t_down - K, 0)
    elif option_type == "put":
        payoff_t_down = np.maximum(K - ST_t_down, 0)
    
    theta = (np.mean(payoff_t_down) - np.mean(payoffs)) / epsilon_t * np.exp(-r * T)
    
    # Vega (sensitivity to volatility)
    epsilon_sigma = 1e-4
    sigma_up = sigma + epsilon_sigma
    ST_sigma_up = S0 * np.exp((r - 0.5 * sigma_up**2) * T + sigma_up * np.sqrt(T) * Z)
    if option_type == "call":
        payoff_sigma_up = np.maximum(ST_sigma_up - K, 0)
    elif option_type == "put":
        payoff_sigma_up = np.maximum(K - ST_sigma_up, 0)
    
    vega = (np.mean(payoff_sigma_up) - np.mean(payoffs)) / epsilon_sigma * np.exp(-r * T)
    
    # Rho (sensitivity to interest rate)
    epsilon_r = 1e-4
    r_up = r + epsilon_r
    ST_r_up = S0 * np.exp((r_up - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    if option_type == "call":
        payoff_r_up = np.maximum(ST_r_up - K, 0)
    elif option_type == "put":
        payoff_r_up = np.maximum(K - ST_r_up, 0)
    
    rho = (np.mean(payoff_r_up) - np.mean(payoffs)) / epsilon_r * np.exp(-r * T)
    
    return option_price, delta, gamma, theta, vega, rho

# Parameters
S0 = 100     # Initial stock price
K = 105      # Strike price
r = 0.05     # Risk-free rate
T = 1        # Time to maturity (1 year)
sigma = 0.2  # Volatility (20%)
num_simulations = 100000  # Number of Monte Carlo simulations

# Call Option Pricing with Full Greeks
call_price, call_delta, call_gamma, call_theta, call_vega, call_rho = monte_carlo_option_pricing_with_full_greeks(
    S0, K, r, T, sigma, num_simulations, "call"
)
print(f"Call Option Price: {call_price:.2f}")
print(f"Delta: {call_delta:.4f}")
print(f"Gamma: {call_gamma:.4f}")
print(f"Theta: {call_theta:.4f}")
print(f"Vega: {call_vega:.4f}")
print(f"Rho: {call_rho:.4f}")

# Put Option Pricing with Full Greeks
put_price, put_delta, put_gamma, put_theta, put_vega, put_rho = monte_carlo_option_pricing_with_full_greeks(
    S0, K, r, T, sigma, num_simulations, "put"
)
print(f"\nPut Option Price: {put_price:.2f}")
print(f"Delta: {put_delta:.4f}")
print(f"Gamma: {put_gamma:.4f}")
print(f"Theta: {put_theta:.4f}")
print(f"Vega: {put_vega:.4f}")
print(f"Rho: {put_rho:.4f}")
