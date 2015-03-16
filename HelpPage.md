NOMBRE
> Caipyranha - Procesador y analizador de secuencias

MODO DE EMPLEO
> python caipyranha.py [-e] [-m] [-a] [-f] [-b]

> Los nombres de los ficheros y carpetas deben ir sin espacios.

> Se debe respetar el orden marcado arriba excepto para la opción -b, pudiéndose saltar alguna
> opción excepto -e.

> Debe indicar la pareja de primers utilizada para amplificar el producto.

> Debe indicar el intervalo de la etiqueta que quiere utilizar para hacer comparaciones.

OPCIONES
> -e 	 Extrae el inserto del plásmido si la puntuación alineamiento local de uno de los primers
> > con la secuencia es mayor del 90% a la máxima. Utiliza la función "pairwise2" de Biopython.

> -m 	 Si dos etiquetas son iguales en el intervalo seleccionado, se unen las secuencia extraídas
> > por las regiones solapantes. Utiliza "merger" de la suite EMBOSS.

> -a 	 Alineamiento de las secuencias extraídas y/o unidas mediante MAFFT.
> -f 	 Obtener árboles filogenéticos mediante los métodos de Neighbor-Joining y Máxima Verosimilitud
> > mediante PHYLIPnew

> -b 	 Realizar un BLAST a través de la web y devuelve el número de resultados indicado.