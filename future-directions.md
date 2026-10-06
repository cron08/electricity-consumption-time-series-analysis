## Future Directions

The current analysis focuses primarily on identifying the temporal and seasonal
structure of electricity consumption, comparing candidate ARIMA/SARIMA models,
and evaluating model adequacy through parameter significance, information
criteria and residual diagnostics. It does not evaluate the models on
out-of-sample forecasting performance.

A natural extension would therefore be to introduce a rolling-origin or
walk-forward forecasting framework and compare models using metrics such as
MAE, RMSE and MASE, alongside prediction-interval coverage. This would test
whether the model preferred on in-sample statistical criteria also provides
the strongest forecasts on unseen observations.

Several additional time-series models could also be investigated:

- **ARFIMA / fractional differencing** to test whether the series exhibits
  persistent long-memory behaviour that is not fully captured by conventional
  integer differencing.
- **Dynamic harmonic regression with Fourier terms and ARIMA errors** to model
  annual seasonality more flexibly than a fixed seasonal ARIMA structure.
- **Exponential smoothing / ETS and state-space models** as alternative
  approaches to modelling trend and seasonal dynamics.
- **TBATS** for flexible modelling of potentially complex or evolving
  seasonality.
- **ARCH/GARCH-type models** if residual diagnostics indicate time-varying
  conditional variance.
- With substantially more data, machine-learning or deep-learning forecasting
  methods could also be benchmarked against the statistical models.

These extensions would shift the project from primarily structural time-series
analysis toward a direct comparison of both explanatory adequacy and
forecasting performance.
