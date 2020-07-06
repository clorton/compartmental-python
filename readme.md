# SSA/Compartmental Modeling with Python

Compare GillesPy2 and StochPy with IDM CMS for disease modeling.

## GillesPy2

[documentation](https://gillespy2.github.io/GillesPy2/docs/build/html/index.html)

## StochPy

[documentation](http://stochpy.sourceforge.net/html/userguide.html)

## IDM CMS

[documentation](https://idmod.org/docs/cms/index.html)  
[repository](https://github.com/institutefordiseasemodeling/idm-cms)

## Test Scenarios

1. simple SEIR model, COVID-19 parameters
2. SEIR model, COVID-19 parameters, with scheduled vaccination
3. SEIR model, COVID-19 parameters, with social distancing triggered by new case counts

|   Species   |   Description            | Initial Value |
|:-----------:|:------------------------:|:-------------:|
| **S**       | **S**usceptible          | 9999          |
| **E**       | **E**xposed/incubating   | 0             |
| **Y**       | s**Y**mptomatic (total)  | 0             |
| **A**       | **A**symptomatic (total) | 0             |
| **I**       | **I**nfectious           | 1             |
| **C**       | **C**ases (total)        | 1             |
| **R**       | **R**ecovered/removed    | 0             |


|  Parameter  | Description                   | Value       |
|:-----------:|-------------------------------|-------------|
|K<sub>i</sub>|infection parameter            | 2.4/6 = 0.4 |
|K<sub>s</sub>|proportion of symptomatic cases| 20% = 0.2   |
|K<sub>e</sub>|incubation parameter           | 1/4 = 0.25  |
|K<sub>r</sub>|recovery parameters            | 1/6 = 0.167 |

```
S → E : Ki * S * I / (S + E + I + R)  
E → I + Y + C : Ks * Ke * E  
E → I + A + C : (1 - Ks) * Ke * E  
I → R : Kr * I
```
