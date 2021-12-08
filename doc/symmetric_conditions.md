Considering symetric boundary conditions, here are all the possible cases: 

![test-2](/uploads/3cdcfee9c2e981ba77e1585a6fa2f622/test-2.png)

This will impact the parameter in the FFTW_Plan_Dim class. 

## Cell centered 
|parameter \ BCs | ODD - ODD              | ODD - EVEN | EVEN - ODD | EVEN - EVEN       |
| ------         | ------                 | ------     | ------     | ------            |
|Transform type  | RODFT10 & RODFT01      | RODFT11    | REDFT11    | REDFT10 & REDFT01 |
|  | | | | |

## Vertex centered 
|parameter \ BCs | ODD - ODD              | ODD - EVEN          | EVEN - ODD            | EVEN - EVEN       |
| ------         | ------                 | ------              | ------                | ------            |
|Transform type  | RODFT00                | RODFT01 & RODFT10   | REDFT01 & REDFT10     | REDFT00           |
|  | | | | |
