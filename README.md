# BAR: nanox maleable

Este repositorio contiene una versión maleable de [nanox](https://github.com/bsc-pm/nanox), que permite modificar el número de hilos y la frecuencia de cada core de forma dinámica.

Para la toma de decisiones, utiliza el co-planificador ([BACO](url-missing)) para la comunicación del estado y la recepción de la lógica de control.

La implementación principal del gestor se encuentra en forma de plugins de planificación dentro de: src/plugins/sched/*


