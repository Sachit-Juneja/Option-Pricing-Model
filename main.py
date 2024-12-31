import streamlit as st
import numpy as np
from backend import european_option_pricing, asian_option_pricing, barrier_option_pricing
from utils import calculate_greeks
from option_visualization import plot_price_paths, plot_payoff, plot_greeks_vs_parameter

st.title("Option Pricing and Greeks Analysis")

# Sidebar for option type selection
option_type = st.sidebar.selectbox("Choose Option Type", ["European", "Asian", "Barrier"])
option_subtype = st.sidebar.selectbox("Call or Put Option", ["call", "put"])

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
        delta, gamma, theta, vega, rho = calculate_greeks(
            european_option_pricing, S0, K, r, T, sigma, num_simulations, option_subtype
        )
        st.success(f"European {option_subtype.capitalize()} Option Price: {price:.2f}")
        st.write(f"Delta: {delta:.4f}, Gamma: {gamma:.4f}, Theta: {theta:.4f}, Vega: {vega:.4f}, Rho: {rho:.4f}")

elif option_type == "Asian":
    if st.button("Calculate Asian Option Price and Greeks"):
        price = asian_option_pricing(S0, K, r, T, sigma, num_simulations, option_subtype)
        delta, gamma, theta, vega, rho = calculate_greeks(
            asian_option_pricing, S0, K, r, T, sigma, num_simulations, option_subtype
        )
        st.success(f"Asian {option_subtype.capitalize()} Option Price: {price:.2f}")
        st.write(f"Delta: {delta:.4f}, Gamma: {gamma:.4f}, Theta: {theta:.4f}, Vega: {vega:.4f}, Rho: {rho:.4f}")

elif option_type == "Barrier":
    barrier = st.number_input("Barrier Level", value=110.0)
    barrier_type = st.selectbox("Barrier Type", ["up-and-out", "down-and-out", "up-and-in", "down-and-in"])
    if st.button("Calculate Barrier Option Price and Greeks"):
        price = barrier_option_pricing(S0, K, r, T, sigma, num_simulations, option_subtype, barrier, barrier_type)
        delta, gamma, theta, vega, rho = calculate_greeks(
            barrier_option_pricing, S0, K, r, T, sigma, num_simulations, option_subtype, barrier=barrier, barrier_type=barrier_type
        )
        st.success(f"Barrier {option_subtype.capitalize()} Option Price ({barrier_type}): {price:.2f}")
        st.write(f"Delta: {delta:.4f}, Gamma: {gamma:.4f}, Theta: {theta:.4f}, Vega: {vega:.4f}, Rho: {rho:.4f}")

# Visualizations
st.sidebar.subheader("Visualization Options")
if st.sidebar.checkbox("Show Simulated Price Paths"):
    num_paths = st.sidebar.slider("Number of Paths to Simulate", min_value=10, max_value=1000, value=100)
    plot_price_paths(S0, r, T, sigma, num_paths)

if st.sidebar.checkbox("Show Option Payoff"):
    plot_payoff(S0, K, option_subtype)

if st.sidebar.checkbox("Show Greek Sensitivity"):
    greek_name = st.sidebar.selectbox("Select Greek", ["Delta", "Gamma", "Theta", "Vega", "Rho"])
    parameter_name = st.sidebar.selectbox("Parameter to Vary", ["S0", "sigma", "r", "T", "K"])
    parameter_min = st.sidebar.number_input(f"Min {parameter_name.capitalize()}", value=0.5 * eval(parameter_name))
    parameter_max = st.sidebar.number_input(f"Max {parameter_name.capitalize()}", value=1.5 * eval(parameter_name))
    parameter_steps = st.sidebar.slider("Steps", min_value=10, max_value=100, value=50)
    parameter_range = np.linspace(parameter_min, parameter_max, parameter_steps)

    # Map Greek names to specific index in calculate_greeks output
    greek_index_map = {"Delta": 0, "Gamma": 1, "Theta": 2, "Vega": 3, "Rho": 4}
    greek_index = greek_index_map[greek_name]

    # Define a function to compute the Greek for a given parameter value
    def greek_function(S0, K, r, T, sigma, parameter_value, parameter_name):
        params = {"S0": S0, "K": K, "r": r, "T": T, "sigma": sigma}
        params[parameter_name] = parameter_value  # Dynamically update the parameter
        return calculate_greeks(
            european_option_pricing,  # Adjust for Asian or Barrier as needed
            params["S0"], params["K"], params["r"], params["T"], params["sigma"],
            num_simulations, option_subtype
        )[greek_index]

    # Pass the function and parameter range to plot_greeks_vs_parameter
    plot_greeks_vs_parameter(S0, K, r, T, sigma, num_simulations, option_subtype,
                             greek_function, parameter_name, parameter_range)
