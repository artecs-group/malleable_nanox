# BAR: nanox maleable

Este repositorio contiene una versión maleable de [nanox](https://github.com/bsc-pm/nanox), que permite modificar el número de hilos y la frecuencia de cada core de forma dinámica.

Para la toma de decisiones, utiliza el co-planificador ([BACO](url-missing)) para la comunicación del estado y la recepción de la lógica de control.

La implementación principal del gestor se encuentra en forma de plugins de planificación dentro de: src/plugins/sched/*


### Acknowledgements
This work was partially supported by the EU (FEDER) and Spanish MINECO (RTI2018-093684-B-I00), MICINN (PID2021-126576NB-I00), CM (S2018/TCS-4423), and the Madrid Government (Comunidad de Madrid- Spain) under the Multiannual Agreement with Complutense University in the line Program to Stimulate Research for Young Doctors in the context of the V PRICIT (Regional Programme of Research and Technological Innovation) under project PR65/19-22445
