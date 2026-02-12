

## Purpose of Volcano-Meteo-Tsurrogate v1.0

The **Volcano-Meteo-Tsurrogate v1.0 software** has been developed to provide **faster-than-real-time forecasts of meteotsunami waves** following explosive volcanic eruptions. Because **large eruptions (VEI ≥ 5)** are rare, the surrogate models implemented in this framework are **trained entirely on deterministic numerical simulations**. To overcome the high computational cost of coupled atmosphere–ocean models, the software employs a **stochastic surrogate modeling strategy based on generalized Polynomial Chaos Expansion (gPCE)**.

This approach enables **efficient uncertainty propagation**, as statistical moments and sensitivity indices can be directly derived from the gPCE coefficients. The **validity** of using surrogate models in this context relies on the **convergence of the gPCE representation** for key hazard metrics (**maximum elevation, maximum current speed, and arrival time**). This assumption has previously been **validated for meteotsunami hazards driven by atmospheric disturbances in the Adriatic Sea** and forms the basis for its **extension to volcano-driven planetary meteotsunamis**.
