# Non-autonomous-problems Solver (Python)

This repository contains standalone Python scripts that solve various non-autonomous differential equations using method of seperation of variable.

## 📂 Examples

- **Ricatti Problem** (`riccati.py`)  
  Nonlinear first-order differential equation, often used in control theory and dynamic systems.
- **Equation 1/(10-t)** (`singularEquation.py`)  
  Singular behavior near t=10; models time-dependent divergence in simple systems.
- **Airy Equation** (`airy.py`)  
  Linear second-order differential equation with applications in quantum mechanics and wave propagation.


## 🧠 Methodology

Each example applies a technique that reduces N-dimensional differential systems into 1D separable parts, then solves them using exponential product formulas (Lie-Trotter, Strang, or Higher order splitting).  


## 📦 Requirements

Install the required packages using:

```bash
pip install -r requirements.txt
```

Contents of reeuirement.txt:
numpy, matplotlib, scipy

## ▶️ How to Run
Open one of the .py files. Inside each file uder 'Initial Parameters', choose your initial values and method of choice (listed). 
```bash
python riccati.py

pyhom singularEquation.py

python airy.py

```

| Model                        | Equation Form                                                                                                                                  | Parameters/Description                              |
| ---------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------- | --------------------------------------------------- |
| **Riccati Equation**      | $\frac{dx}{dt} = x^2 + t^2$, $x(s)=x_0$ | Nonlinear first-order ODE; often used in control theory and dynamic systems |
| **Singular Equation**           | $\frac{dx}{dt} = \frac{1}{10 - t}$, $x(s)=x_0$ |Exhibits singularity at `t = 10`; useful for testing numerical stability near divergence |
| **Airy Equation**    |  $\frac{d^2y}{dt^2} - t y = 0$, $y(s)=x$, $y'(s)=z$ |  Linear second-order ODE; arises in quantum mechanics and wave propagation problems    |



## 📊 Output

Each script may output:

  Solution curves as plots (saved as PNGs or JPGs in results/)

  Comparison graph with RK45 

  You can modify initial values to save or display these results as needed.

## 📌 Comparing RK45  

We implement the Runge-Kutta RK45 method for solving ODEs. You can find a detailed breakdown of the equations, step-by-step methodology, and an example problem in [Compare RK45](Compare_RK45.md).


## 📌 License

This project is licensed under the [MIT License](LICENSE). You’re free to use, modify, and distribute it with attribution.

## 🙋‍♂️ Authors

**Arun Banjara**  
Ph.D. in Mathematics | Louisiana State University (LSU) 

**Ibrahem ALJabea**  
*Co-author*  
Ph.D. in Mathematics | Louisiana State University (LSU)

## 🔗 Related Projects
If you’d like to see a web version of these solvers with interactivity and graph comparisons, visit the Flask Web App version
[ODE SOLVER](https://arun1111.pythonanywhere.com/)
