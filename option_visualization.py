import numpy as np
import matplotlib.pyplot as plt
import streamlit as st

def plot_price_paths(S0, r, T, sigma, num_paths):
    dt = T / 100  # Time step
    paths = np.zeros((101, num_paths))
    paths[0] = S0
    for t in range(1, 101):
        Z = np.random.standard_normal(num_paths)
        paths[t] = paths[t - 1] * np.exp((r - 0.5 * sigma**2) * dt + sigma * np.sqrt(dt) * Z)
    
    plt.figure(figsize=(10, 6))
    for i in range(min(num_paths, 10)):  # Plot up to 10 paths for clarity
        plt.plot(np.linspace(0, T, 101), paths[:, i])
    plt.title("Simulated Stock Price Paths")
    plt.xlabel("Time to Maturity")
    plt.ylabel("Stock Price")
    plt.grid()
    st.pyplot(plt)

def plot_payoff(S0, K, option_type):
    ST = np.linspace(0.5 * S0, 1.5 * S0, 100)
    if option_type == "call":
        payoffs = np.maximum(ST - K, 0)
    else:
        payoffs = np.maximum(K - ST, 0)
    
    plt.figure(figsize=(10, 6))
    plt.plot(ST, payoffs, label=f"{option_type.capitalize()} Option Payoff")
    plt.axvline(x=K, color="red", linestyle="--", label="Strike Price")
    plt.title(f"{option_type.capitalize()} Option Payoff")
    plt.xlabel("Stock Price at Maturity (ST)")
    plt.ylabel("Payoff")
    plt.legend()
    plt.grid()
    st.pyplot(plt)

def plot_greeks_vs_parameter(S0, K, r, T, sigma, num_simulations, option_type, greek_function, parameter_name, parameter_range):
    greek_values = []
    for parameter in parameter_range:
        kwargs = {parameter_name: parameter}
        greek = greek_function(S0, K, r, T, sigma, num_simulations, option_type, **kwargs)
        greek_values.append(greek)
    
    plt.figure(figsize=(10, 6))
    plt.plot(parameter_range, greek_values, label=f"{parameter_name.capitalize()} Sensitivity")
    plt.title(f"Greek Sensitivity to {parameter_name.capitalize()}")
    plt.xlabel(parameter_name.capitalize())
    plt.ylabel("Greek Value")
    plt.legend()
    plt.grid()
    st.pyplot(plt)
