import numpy as np

import numpy as np

def monte_carlo_option_pricing(S0, K, r, T, sigma, num_simulations):
    """
    Monte Carlo simulation for European call option pricing.

    Parameters:
    - S0: Initial stock price
    - K: Strike price
    - r: Risk-free interest rate
    - T: Time to maturity (in years)
    - sigma: Volatility of the underlying asset
    - num_simulations: Number of simulations

    Returns:
    - Estimated price of the call option
    """
    # Generate random standard normal numbers
    Z = np.random.standard_normal(num_simulations)
    
    # Simulate end prices for the stock
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    
    # Calculate the payoff for a call option
    payoffs = np.maximum(ST - K, 0)
    
    # Discount payoffs back to present value
    discounted_payoffs = payoffs * np.exp(-r * T)
    
    # Average the discounted payoffs to estimate the option price
    option_price = np.mean(discounted_payoffs)
    
    return option_price

# Parameters
S0 = 100     # Initial stock price
K = 105      # Strike price
r = 0.05     # Risk-free rate
T = 1        # Time to maturity (1 year)
sigma = 0.2  # Volatility (20%)
num_simulations = 100000  # Number of Monte Carlo simulations

# Estimate the option price
option_price = monte_carlo_option_pricing(S0, K, r, T, sigma, num_simulations)
print(f"Estimated European Call Option Price: {option_price:.2f}")
