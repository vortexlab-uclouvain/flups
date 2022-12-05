Considering symetric boundary conditions, here are all the possible cases: 

![test-2](/uploads/3cdcfee9c2e981ba77e1585a6fa2f622/test-2.png)

This will impact the parameter in the FFTW_Plan_Dim class. 

## Cell centered RHS 
|parameter \ BCs       | ODD - ODD              | ODD - EVEN        | EVEN - ODD        |  EVEN - EVEN          |
| ------               | :------:               | :------:          | :------:          | :------:              |
| Transform type       | RODFT10 & RODFT01      | RODFT11           | REDFT11           | REDFT10 & REDFT01     |
| Normalisation factor | $\frac{1}{2N}$         | $\frac{1}{2N}$    | $\frac{1}{2N}$    | $\frac{1}{2N}$        |
| Offset of the modes  | 0.0                    | 0.5               | 0.5               | 0.0                   |
| # of points out      | N + 1                  | N                 | N                 | N + 1                 |
| Output correction    | The 0 mode = 0         | /                 | /                 | The fli-flop mode = 0 |

## Vertex centered 
|parameter \ BCs       | ODD - ODD                   | ODD - EVEN        | EVEN - ODD        | EVEN - EVEN           |
| ------               | :------:                    | :------:          | :------:          | :------:              |
|Transform type        | RODFT00                     | RODFT01 & RODFT10 | REDFT01 & REDFT10 | REDFT00               |
| Normalisation factor | $\frac{1}{2(N+1)}$          | $\frac{1}{2N}$    | $\frac{1}{2N}$    | $\frac{1}{2(N-1)}$    |
| Offset of the modes  | 0.0                         | 0.0               | 0.0               | 0.0                   |
| # of points out      | N - 2                       | N - 1             | N - 1             | N                     |
| Input correction     | First and last data ignored | First data ignored| Last data ignored | /                     |