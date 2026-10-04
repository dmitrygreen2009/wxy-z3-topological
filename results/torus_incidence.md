# Four-star periodic quotient

Primitive honeycomb Bravais coordinates `(x,y)`, with A joined to B at
`(x,y)`, `(x-1,y)`, `(x,y-1)` for legs 1, 2, 3. Periods `(2,0)`, `(0,1)`.
Vertex order: A(0,0), B(0,0), A(1,0), B(1,0).

| Physical edge (one-based) | Endpoint 1 | Endpoint 2 | Local leg at both endpoints |
|---|---|---|---|
|1|A(0,0)|B(0,0)|1|
|2|A(0,0)|B(1,0)|2|
|3|A(0,0)|B(0,0)|3|
|4|A(1,0)|B(1,0)|1|
|5|A(1,0)|B(0,0)|2|
|6|A(1,0)|B(1,0)|3|

Local-leg table: `[1,2,3]`, `[1,5,3]`, `[4,5,6]`, `[4,2,6]`.
Matter site indices are 1–12, gauge site indices are 13–18.
Parallel edges 1/3 and 4/6 are distinct physical spins. This multigraph is not K4.
