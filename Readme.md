Projektet ''Indsivning af uvedkommende vand i afløbssystemer'' er et afgangsspeciale på 3. - 4. semester af kandidaten i 'Water and Environmental Engineering' (Vand og Miljø) på Aalborg Universitet.

Den numeriske model, som sammenligner med en eksperimentel sandboks, regner på kvantificering og tidspunkt for indsivning i et rørsystem.
Modellen regner i 1D, med et rørindsivende pseudo 2D-element og tager derfor ikke højde for sandboksens horisontale strømninger. Modellen rammer tidspunktet for indsivning tilsvarende med hvad der blev observeret,
men på grund af 1D-modellering, indsiver fem gange så meget som det målte.

Modellen kan regne vertikalt ned, horisontalt eller vertikalt op, hvor den i rapporten har fokuseret på de vertikale strømninger.
Modellen loader en hotstart ('Hot_start.npy') som initialbetingelse, hvor modellen regner vertikalt ned med denne hotstart.
Hotstarten er resultatet af en 72 timers simulering i den vertikale retning mod gravitation, for at indstille jordsøjlen i en ligevægt og opnå et kapillært vandspejl.

Efter en simulering gemmes NumPy filer, der senere bruges til at visualisere resultaterne.

