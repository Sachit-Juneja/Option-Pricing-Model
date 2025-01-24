# Option Pricing and Risk Analysis with ML and AI

## Overview
This project is a comprehensive tool for option pricing, risk analysis, and financial modeling. It incorporates traditional mathematical models, such as Black-Scholes and Monte Carlo simulations, alongside advanced machine learning (ML) and artificial intelligence (AI) techniques to analyze and predict option prices and sensitivities (Greeks). The tool also includes interactive visualizations and an intuitive interface for exploring complex financial instruments.

## Features
1. **Traditional Option Pricing Models**
   - Black-Scholes for European options.
   - Monte Carlo simulations for Asian and Barrier options.
2. **Advanced Greeks Calculation**
   - Delta, Gamma, Vega, Theta, and Rho for different option types.
   - Sensitivity analysis for portfolio risk management.
3. **Machine Learning Models**
   - Predict option prices using regression models and neural networks.
   - Feature engineering for volatility, time-to-maturity, and more.
4. **Interactive Visualizations**
   - Payoff diagrams and sensitivity heatmaps.
   - Greek sensitivity graphs and implied volatility surfaces.
5. **User-Friendly Front-End**
   - Built with Streamlit for interactive user input and dynamic results.


## Installation
1. Clone this repository:
   ```bash
   git clone https://github.com/your-username/option-pricing-modell.git
   ```
2. Navigate to the project directory:
   ```bash
   cd option-pricing-ml
   ```
3. Create a virtual environment:
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows, use `venv\Scripts\activate`
   ```
4. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage
Run the Streamlit app:
```bash
streamlit run src/frontend/streamlit_app.py
```
This will open an interactive web app in your browser where you can:
- Choose an option type (European, Asian, Barrier).
- Input parameters such as strike price, volatility, and time-to-maturity.
- View calculated prices, Greeks, and sensitivity visualizations.

## Key Technologies
- **Programming Languages**: Python
- **Libraries and Frameworks**:
  - Numerical Computations: `numpy`, `scipy`
  - Machine Learning: `scikit-learn`, `tensorflow`, `xgboost`
  - Visualization: `matplotlib`, `seaborn`, `plotly`
  - Web Interface: `streamlit`

## Features to Explore
1. **Payoff Diagrams**
   - Visualize the profit and loss (P&L) for European, Asian, and Barrier options.
2. **Sensitivity Graphs**
   - Observe how Greeks vary with changes in parameters like underlying price and volatility.
3. **Machine Learning Models**
   - Compare predicted option prices with traditional methods.
4. **Portfolio Analysis**
   - Aggregate Greeks for multiple options to analyze risk.

## Future Enhancements
- Integration of real-time market data using APIs (e.g., Yahoo Finance, Alpha Vantage).
- Reinforcement learning for dynamic hedging strategies.
- Calibration of advanced stochastic models (e.g., Heston, SABR).
- Backtesting framework for evaluating option trading strategies.

## Contributing
Contributions are welcome! Feel free to open an issue or submit a pull request if you have ideas to improve this project.


