import streamlit as st
import numpy as np

def european_option_pricing(S0, K, r, T, sigma, num_simulations, option_type):
    Z = np.random.standard_normal(num_simulations)
    ST = S0 * np.exp((r - 0.5 * sigma**2) * T + sigma * np.sqrt(T) * Z)
    if option_type == "call":
        payoffs = np.maximum(ST - K, 0)
    else:
        payoffs = np.maximum(K - ST, 0)
    price = np.mean(payoffs) * np.exp(-r * T)
    return price

def asian_option_pricing(S0, K, r, T, sigma, num_simulations, option_type):
    dt = T / 100  # Divide the time to maturity into 100 steps
    ST_path = S0 * np.cumprod(np.exp((r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * np.random.randn(100, num_simulations)), axis=0)
    avg_price = np.mean(ST_path, axis=0)
    if option_type == "call":
        payoffs = np.maximum(avg_price - K, 0)
    else:
        payoffs = np.maximum(K - avg_price, 0)
    price = np.mean(payoffs) * np.exp(-r * T)
    return price

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
    price = np.mean(payoffs) * np.exp(-r * T)
    return price

def calculate_greeks(pricing_function, S0, K, r, T, sigma, num_simulations, option_type, **kwargs):
    epsilon = 1e-4  # Small change for numerical derivatives
    
    # Price at S0
    price = pricing_function(S0, K, r, T, sigma, num_simulations, option_type, **kwargs)
    
    # Delta: Sensitivity to S0
    price_up = pricing_function(S0 + epsilon, K, r, T, sigma, num_simulations, option_type, **kwargs)
    price_down = pricing_function(S0 - epsilon, K, r, T, sigma, num_simulations, option_type, **kwargs)
    delta = (price_up - price_down) / (2 * epsilon)
    
    # Gamma: Sensitivity of Delta to S0
    gamma = (price_up - 2 * price + price_down) / (epsilon ** 2)
    
    # Theta: Sensitivity to T (time decay)
    price_shorter = pricing_function(S0, K, r, T - epsilon, sigma, num_simulations, option_type, **kwargs)
    theta = (price_shorter - price) / epsilon
    
    # Vega: Sensitivity to σ (volatility)
    price_vol_up = pricing_function(S0, K, r, T, sigma + epsilon, num_simulations, option_type, **kwargs)
    vega = (price_vol_up - price) / epsilon
    
    # Rho: Sensitivity to r (risk-free rate)
    price_rate_up = pricing_function(S0, K, r + epsilon, T, sigma, num_simulations, option_type, **kwargs)
    rho = (price_rate_up - price) / epsilon
    
    return delta, gamma, theta, vega, rho

# Streamlit UI
st.title("Option Pricing with Greeks")

option_type = st.sidebar.selectbox("Choose Option Type", ["European", "Asian", "Barrier"])
option_subtype = st.sidebar.selectbox("Choose Call or Put", ["call", "put"])

# Common inputs
S0 = st.number_input("Initial Stock Price (S0)", value=100.0)
K = st.number_input("Strike Price (K)", value=100.0)
r = st.number_input("Risk-Free Rate (r, in %)", value=5.0) / 100
T = st.number_input("Time to Maturity (T, in years)", value=1.0)
sigma = st.number_input("Volatility (σ, in %)", value=20.0) / 100
num_simulations = st.number_input("Number of Simulations", value=100000, step=1000)

if option_type == "European":
    if st.button("Calculate European Option Price and Greeks"):
        price = european_option_pricing(S0, K, r, T, sigma, num_simulations, option_subtype)
        delta, gamma, theta, vega, rho = calculate_greeks(european_option_pricing, S0, K, r, T, sigma, num_simulations, option_subtype)
        
        st.success(f"European {option_subtype.capitalize()} Option Price: {price:.2f}")
        st.write(f"**Delta**: {delta:.4f}")
        st.write(f"**Gamma**: {gamma:.4f}")
        st.write(f"**Theta**: {theta:.4f}")
        st.write(f"**Vega**: {vega:.4f}")
        st.write(f"**Rho**: {rho:.4f}")

elif option_type == "Asian":
    if st.button("Calculate Asian Option Price and Greeks"):
        price = asian_option_pricing(S0, K, r, T, sigma, num_simulations, option_subtype)
        delta, gamma, theta, vega, rho = calculate_greeks(asian_option_pricing, S0, K, r, T, sigma, num_simulations, option_subtype)
        
        st.success(f"Asian {option_subtype.capitalize()} Option Price: {price:.2f}")
        st.write(f"**Delta**: {delta:.4f}")
        st.write(f"**Gamma**: {gamma:.4f}")
        st.write(f"**Theta**: {theta:.4f}")
        st.write(f"**Vega**: {vega:.4f}")
        st.write(f"**Rho**: {rho:.4f}")

elif option_type == "Barrier":
    barrier = st.number_input("Barrier Level", value=110.0)
    barrier_type = st.selectbox("Barrier Type", ["up-and-out", "down-and-out", "up-and-in", "down-and-in"])
    if st.button("Calculate Barrier Option Price and Greeks"):
        price = barrier_option_pricing(S0, K, r, T, sigma, num_simulations, option_subtype, barrier, barrier_type)
        delta, gamma, theta, vega, rho = calculate_greeks(barrier_option_pricing, S0, K, r, T, sigma, num_simulations, option_subtype, barrier=barrier, barrier_type=barrier_type)
        
        st.success(f"Barrier {option_subtype.capitalize()} Option Price ({barrier_type}): {price:.2f}")
        st.write(f"**Delta**: {delta:.4f}")
        st.write(f"**Gamma**: {gamma:.4f}")
        st.write(f"**Theta**: {theta:.4f}")
        st.write(f"**Vega**: {vega:.4f}")
        st.write(f"**Rho**: {rho:.4f}")
