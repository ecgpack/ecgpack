# Input File Format

Table of Contents:

- [Basic structure of the input file](#basic-structure-of-the-input-file)
- [Input file sections](#input-file-sections)
  - [Header](#header)
  - [Command list](#command-list)
  - [History](#history)
  - [Basis functions](#basis-functions)
- [Creating input files](#creating-input-files)
- [Restarting calculations that generate or optimize a basis](#restarting-calculations-that-generate-or-optimize-a-basis)

## A short example of an input file

An example of an input file `inout.txt` for the ground state calculation of He is provided below. This input file contains only 10 basis functions.

```
 BASIS_TYPE RG_0S
 PARTICLES      3
 MASSES  0.7294299541700000E+04  0.1000000000000000E+01  0.1000000000000000E+01
 CHARGES  0.2000000000000000E+01 -0.1000000000000000E+01 -0.1000000000000000E+01
 SYMMETRY (1+P23)
 BASIS_SIZE     10
 CURRENT_ENERGY -0.2900744334456869E+01
 WHICH_EIGENVALUE      1
 EIGVAL_TOLERANCE  0.1000000000000000E-11
 INVITPARAMETER  0.1000001000000000E+01
 LAST_EIGVAL_TOL  0.1751847624688021E-13
 BEST_EIGVAL_TOL  0.4290247453116364E-14
 WORST_EIGVAL_TOL  0.9986026673536815E-12
 GENERATOR_PARAM  0.7000000000000000E+00  0.1000000000000000E+01  0.3000000000000000E+01
 ==============================
     10 -0.2900744334456869E+01     10      0    123
 ==============================
 BASIS_ENL G      1      5      1    500    200  0.9500000000000000E+00  0.3000000000000000E+01
 FULL_OPT1 G      5      1      5 999999  0.9500000000000000E+00  0.3000000000000000E+01     60     60  none
 FULL_OPT1 G      5      1      5 999999  0.9500000000000000E+00  0.3000000000000000E+01     60     60  none
 OPT_CYCLE G      5      1      5      1      1     10    130  0.9500000000000000E+00  0.3000000000000000E+01     10
 SAVE_FILE      5  inout_He_1Se-01-00005.txt
 BASIS_ENL I      6     10      1    500    200  0.9500000000000000E+00  0.3000000000000000E+01
 FULL_OPT1 I     10      1     10 999999  0.9500000000000000E+00  0.3000000000000000E+01     60     60  none
 OPT_CYCLE I     10      1     10      1      1     10    130  0.9500000000000000E+00  0.3000000000000000E+01     10
 SAVE_FILE     10  inout_He_1Se-01-00010.txt
 ==============================
      1 -0.1668865998401775E+01      0      0      0
      2 -0.2642629552717170E+01      0      0      0
      3 -0.2798648895724608E+01      0      0      0
      4 -0.2844201302182638E+01      0      0      0
      5 -0.2889258918771523E+01     10      0    121
      6 -0.2893514456170832E+01      0      0      0
      7 -0.2895565511375141E+01      0      0      0
      8 -0.2896518341033143E+01      0      0      0
      9 -0.2897179931311348E+01      0      0      0
     10 -0.2900744334456869E+01     10      0    123
 ==============================
      1  0.8225526505330379E+00  0.3286331967568486E-01  0.4734092048850279E+00
      2  0.4657831516551231E+00 -0.4035888660117890E-01  0.1841499686869463E+01
      3  0.5792314542342561E+00 -0.2455636629377815E-01  0.8488492583438731E+01
      4  0.1316320598387222E+01  0.1264307213755476E+00  0.7672315752278411E+00
      5  0.1230152783754290E+01  0.1911586281826508E+00  0.2384612056831933E+01
      6 -0.2233523534916765E+01  0.1516523349832122E+01 -0.1098568271679442E+01
      7  0.9104134991273443E+00  0.2902311989845964E-01  0.1975856897728040E+02
      8  0.1405523968799062E+01  0.3635708756466438E+00  0.6425523571898003E+00
      9 -0.6106458148603032E+00 -0.6976991906481144E-02 -0.3667310281682786E+01
     10  0.1577666796470805E+01  0.2898602080507635E+00  0.5856740074044248E+01
```

## Basic structure of the input file

The input file `inout.txt` consists of four sections: header, list of commands, calculation history, and parameters of basis functions:

| Section | Short Description |
| :--- | :--- |
| **Header** | Defines the quantum system - number of particles, their charges, masses, permutational symmetry. Also contains the current size of the basis, current energy, parameters for the eigensolver and its statistics, as well as parameters of the distribution used in stochastic selection of variational parameters. |
| **Command List** | Contains a list of actions that need to be performed. |
| **History** | Lists the energy values that were obtained at each basis size when the basis was being built, as well as the number of optimization steps that were performed (e.g. number of objective function evaluations, number of optimization cycles). |
| **Basis functions** | Contains parameters of basis functions, one function per line. |

When calculations begin from scratch, only the header and command list need to be supplied in the input file.

Note that in any calculation that involves generation of a basis (either from scratch or by extending an existing basis) the input file also serves as an output file (hence the name - `inout.txt`). This file is updated with any incremental changes of the basis. The frequency with which the file gets overwritten can be controlled by the user.

In many practical situations one must save intermediate basis sets. This can be easily achieved with a command `SAVE_FILE`, which is described below.

## Input file sections

### Header

A typical header looks as follows (it is an example for the ground state calculations of lithium atom that begin from scratch - no basis is yet generated):

```
 BASIS_TYPE RG_0S
 PARTICLES    4
 MASSES    12786.39227775 1.0 1.0 1.0
 CHARGES   3.0 -1.0 -1.0 -1.0
 SYMMETRY  (1+P23)(1-P24)
 BASIS_SIZE  0
 CURRENT_ENERGY    1.0E30
 WHICH_EIGENVALUE  1    
 EIGVAL_TOLERANCE  1.0E-12
 INVITPARAMETER    1.000001000E0
 LAST_EIGVAL_TOL   1.0E0
 BEST_EIGVAL_TOL   1.0E0
 WORST_EIGVAL_TOL  1.0E-30
 GENERATOR_PARAM   0.7E0  1.0E0  3.0E0
 ==============================
    100 -0.7477376516187604E+01     10      0      0
 ==============================
```

At its end, the header includes an excerpt of the last line in the history, surrounded by two ruler lines.
This excerpt shows the current state of the optimization process.

Note that each line of the header (except for the history excerpt) begins with a keyword. Some lines in the header are optional and can be skipped. However, most are required. Lines with the keywords must appear in a specific order and no blank lines are allowed.
The description of all header keywords, in the order they appear in input file, is provided in the table below. Examples of complete input file lines for each keyword are also provided.

| Keyword | Requirement | Description |
| :--- | :---: | :--- |
| `BASIS_TYPE` | optional | Specifies the type of basis. Possible values are `CG_0S`, `RG_0S`, `RG_1P`, `RG_2D`, `RG_2P`. While this line is optional, it is highly recommended to include it in the input file to avoid confusion when running multiple calculations that use different basis types. <br> *Example* : `BASIS_TYPE RG_2D` |
| `PARTICLES` | required | Specifies the number of particles in the system. For example, for Li atom it is 4 (a nucleus + three electrons). <br> *Example* : `PARTICLES 5` |
| `FIXED_INDEX` | optional | This optional line can only be used in the case of `RG_1P` basis type. The keyword is followed by the fixed value of the z-index that should be used in calculations. <br> *Example* : `FIXED_INDEX 2` |
| `FIXED_INDEX_1` | optional | This optional line can only be used in the case of `RG_2D` or `RG_2P` basis types. The keyword is followed by the fixed value of the first integer index that should be used in calculations. <br> *Example* : `FIXED_INDEX_1 2` |
| `FIXED_INDEX_2` | optional | This optional line can only be used in the case of `RG_2D` or `RG_2P` basis types. The keyword is followed by the fixed value of the second integer index that should be used in calculations. <br> *Example* : `FIXED_INDEX_2 3` |
| `MASSES` | required | Specifies masses of all particles (in atomic units, so that mass of the electron is 1.0), starting with the first one. The values are separated by one or more spaces. For calculations of atoms or ions with a clamped nucleus the nuclear mass should be set to a very large value, e.g. 1.0E+30. <br> *Example* : `MASSES 1.0E30 1.0 1.0 1.0 1.0` |
| `CHARGES` | required | Specifies Coulomb charges of all particles (in atomic units), starting with the first one. The values are separated by one or more spaces. It is possible to use noninteger charges. <br> *Example* : `CHARGES 4.0 -1.0 -1.0 -1.0 -1.0` |
| `REPULSION_SCALING_PARAM` | Optional | Scaling parameter for all repulsive interactions in the system (default: `1.0`). When supplied, all positive charge products $q_i q_j$ in both the nonrelativistic and leading relativistic Hamiltonians are scaled as $q_i q_j \rightarrow \gamma q_i q_j$. Setting $\gamma < 1$ is useful for studying the stability of Coulomb systems: one can artificially lower repulsion to establish binding and generate a high-quality wave function, then gradually scale $\gamma$ back to `1.0`. The advantage of scaling all attractive or all repulsive interactions in the system simultaneously in comparison with scaling individual charges or masses is that the former simultaneous scaling does not break permutational or charge conjugation symmetry. <br> *Example* : `REPULSION_SCALING_PARAM 0.99`. |
| `REPULSION_SCALING_PARAM_PLUS` | optional | This scaling parameter (default: `1.0`) is similar to `REPULSION_SCALING_PARAM`, but it scales only the repulsive interactions where both charges $q_i$ and $q_j$ are positive. <br> *Example* : `REPULSION_SCALING_PARAM_PLUS 0.99` |
| `REPULSION_SCALING_PARAM_MINUS` | optional | This scaling parameter (default: `1.0`) is similar to `REPULSION_SCALING_PARAM`, but it scales only the repulsive interactions where both charges $q_i$ and $q_j$ are negative. Note that all three scaling parameters, `REPULSION_SCALING_PARAM`, `REPULSION_SCALING_PARAM_PLUS`, and `REPULSION_SCALING_PARAM_MINUS`, scale charge products independently, so all three can be applied. <br> *Example* : `REPULSION_SCALING_PARAM_MINUS 0.99` |
| `ATTRACTION_SCALING_PARAM` | optional | Scaling parameter $\beta$ for all attractive interactions in the system (default: `1.0`). It works pretty much in the same way as `REPULSION_SCALING_PARAM` but for all negative/attractive products $q_i q_j$ in the system, i.e. it scales $q_i q_j \rightarrow \beta q_i q_j$.  <br> *Example* : `ATTRACTION_SCALING_PARAM 1.01`. |
| `SYMMETRY` | required | The keyword is followed by a string that defines the Young projection operator for systems with identical particles. ECGPACK adopts spin-free formalism in the nonrelativistic calculations of the energy. The Young operator is a product of symmetrizers over all rows and antisymmetrizers over all columns of a suitable Young tableau. For the ground doublet state of lithium, a suitable Young diagram is $[ 2 \ 1]$ , which yields a symmetrizer over the first row `(1+P23)` and an antisymmetrizer over the first column `(1-P24)`. Here `Pij` denotes a pair permutation (transposition) of particles i and j. The Young operator becomes progressively more complicated for systems with a larger number of identical particles, but the structure in the form of the symmetrizer and antisymmetrizer (or, alternatively, in reverse order) remains. For example, for the doublet states of boron (a six-particle system; nucleus is particle 1, electrons are particles 2 through 6) with a Young diagram $[ 2^2 \ 1 ]$  a suitable expression for the Young operator is `(1+P23)(1+P45)(1-P24)(1-P26-P46)(1-P35)`. For the quartet states of boron with a Young diagram $[ 2 \ 1^3]$ the expression for the Young operator is `(1+P23)(1-P24)(1-P25-P45)(1-P26-P46-P56)`. If there is more than one set of identical particles in the system, the total Young operator is a product of Young operators for each set. For some systems there may be additional symmetry - e.g. charge conjugation - that can be enforced by choosing a proper expression for the Young operator. For example, for the ground state of the positronium molecule, where the first two particles are electrons and the third and fourth particles are positrons, a suitable expression is `(1+P13P24)(1+P34)(1+P12)`. <br> *Example* : `SYMMETRY (1+P23)(1-P24)(1-P25-P45)` |
| `BASIS_SIZE` | required | Indicates the current size of the basis. In new calculations it should be set to zero. <br> *Example* : `BASIS_SIZE  100`. |
| `CURRENT_ENERGY` | required | The value of the current energy for the basis stored in the input file. In new calculations that start from scratch this should be set to some very large value, e.g. 1.0E30. <br> *Example* : `CURRENT_ENERGY -0.7477376516187604E+01` |
| `WHICH_EIGENVALUE` | required | Specifies which energy eigenvalue should be targeted in the calculations. 1 means the lowest, 2 means the second lowest, etc. Note that this parameter matters only if the LAPACK eigenvalue solver is used (see the next subsection for more details on that). If the inverse iteration eigenvalue solver is used this parameter can be anything. Also note that if the LAPACK eigensolver is used the value specified by `WHICH_EIGENVALUE` cannot exceed the current size of the basis. Therefore, when initiating calculations of the $k$-th excited state one must first generate or copy a basis containing at least $k+1$ ECG basis functions for the ground (or some lower excited state) and then switch to the eigenvalue number $k$ . <br> *Example* : `WHICH_EIGENVALUE 2` |
| `EIGVAL_TOLERANCE` | required | Specifies the desired accuracy of the inverse iteration eigenvalue solver. It is recommended to make `EIGVAL_TOLERANCE` as small as possible in each particular calculation but not too close to the value of the machine epsilon. In practice, due to rounding errors and some degree of linear dependencies in the generated bases, the value of `EIGVAL_TOLERANCE` should be a couple of orders of magnitude larger than the machine epsilon. Otherwise the inverse iteration eigenvalue solver may fail too often. For calculations with double precision a reasonable starting choice for `EIGVAL_TOLERANCE` is somewhere around 1.0E-12 and as the basis generation progresses it may need to be adjusted to balance between accuracy and stability. When the LAPACK eigenvalue solver is invoked in calculations, the value `EIGVAL_TOLERANCE` is not referenced and does not affect anything. <br> *Example* : `EIGVAL_TOLERANCE 1.0E-12` |
| `INVITPARAMETER` | required | The parameter which defines how close the approximate eigenvalue used in the inverse iteration eigenvalue solver should be to the current energy value. This parameter should normally be very close to `1.0` so that the eigensolver targets the right eigenvalue and there is no accidental root switching. Also, the closer it is to `1.0` the fewer iterations are needed to converge the eigenvector/eigenvalue to the desired accuracy (defined by `EIGVAL_TOLERANCE`). Yet, it must not be exactly `1.0` because if the approximate eigenvalue happens to be equal to the exact eigenvalue then it may lead to a failure of the inverse iteration algorithm. Typically, a value of `1.000001` or `1.0000001` provides a good balance between speed and stability. However, in some calculations it may need to be adjusted. It should be noted that the $H - \varepsilon_{appr} S$ matrix factorization (which involves $\varepsilon_{appr}$ because it depends on `INVITPARAMETER`) is changed every time a new command from the command list (see next subsection) is executed. <br> *Example* : `INVITPARAMETER 1.000001` |
| `LAST_EIGVAL_TOL` | required | The value of the actual relative accuracy of the eigenvalue obtained in the last call of the inverse iteration eigenvalue solver. This value does not change anything and is used for monitoring only. <br> *Example* : `LAST_EIGVAL_TOL 0.1170909566892634E-13` |
| `BEST_EIGVAL_TOL` | required | The value of the best relative accuracy obtained by the inverse iteration eigenvalue solver in the course of the entire calculation since it started from scratch or since this value was reset. This value does not change anything and is used for monitoring only. <br> *Example* : `BEST_EIGVAL_TOL 0.2351296478907559E-14` |
| `WORST_EIGVAL_TOL` | required | The value of the worst relative accuracy obtained by the inverse iteration eigenvalue solver in the course of the entire calculation since it started from scratch or since this value was reset. This value does not change anything and is used for monitoring only. <br> *Example* : `WORST_EIGVAL_TOL 0.9999931369150333E-12` |
| `GENERATOR_PARAM` | required | List of three parameters (let us call them $p_1$ , $p_2$ , and $p_3$ ) which define the shape of the master distribution used at the stage of stochastic selection of new basis functions. This master distribution is a sum of two distributions: one comes with a weight factor of $p_1$ (so that the first supplied parameter must lie between `0.0` and `1.0`), the other with a weight factor $1-p_1$ . The first distribution samples independently each nonlinear variational parameter $\alpha_i$ of existing basis functions (based on a normal distribution centered at $a_i$ and width given by $a_i p_2$ ). The second distribution samples all nonlinear parameters of a new basis function by scaling the parameters of an existing function by $s$ , where $s$ is distributed normally with a center at $1.0$ and width given by $p_3$ (negative values of $s$ are discarded). <br> *Example* :  `GENERATOR_PARAM  0.7 1.0 3.0` |

### Command list

After the header, an input file contains a list of commands that instruct the program which actions need to be performed. For the energy codes most commonly this list includes commands that expand the basis, do cyclic optimization of nonlinear variational parameters, do full optimization of all nonlinear variational parameters, save the basis into a separate file when a certain basis size is reached, and compute expectation values. Here is an example of the command list taken from the same input file for the He atom ground state (some real numbers have been shortened for easier reading):

```
 BASIS_ENL  G    1    5    1  500  200  0.95  3.0
 FULL_OPT1  G    5    1    5 999999  0.95  3.0     60     60  none
 FULL_OPT1  G    5    1    5 999999  0.95  3.0     60     60  none
 OPT_CYCLE  G    5    1    5    1    1   10  130  0.95  3.0     10
 SAVE_FILE     5  inout_He_1Se-01-00005.txt
 BASIS_ENL  I    6   10    1  500    200  0.95  3.0
 FULL_OPT1  I   10    1   10 999999  0.95  3.0   60   60  none
 OPT_CYCLE  I   10    1   10    1    1   10  130  0.95  3.0     10
 SAVE_FILE    10  inout_He_1Se-01-00010.txt
```

Each command contains a number of integer, real, character, or string arguments that must be arranged in a certain order. Spacings between those arguments can be arbitrary. Detailed descriptions of all available commands and their arguments are provided below.

Most commands begin with a one-character **eigenvalue solver type**, which can be either `G` or `I`:

- `G` selects the LAPACK routine `DSYGVX`, which first reduces the definite generalized symmetric eigenvalue problem to a standard symmetric eigenvalue problem and then solves it. When `G` is used, the targeted root is set by the header keyword `WHICH_EIGENVALUE` (the header values `CURRENT_ENERGY` and `EIGVAL_TOLERANCE` are not referenced). The `G` option for the eigenvalue solver is safe in the sense that it prevents unintended root switching (it always targets a specific eigenvalue), which is particularly important for the case of small basis sizes, when the total wave function may still undergo considerable change in the process of its optimization. However, the `G` option is extremely slow when it comes to updating the solution in the case of large basis sizes, which is what is done routinely in `BASIS_ENL`, `OPT_CYCLE`, and `FULL_OPT1` (see the descriptions of these commands below). Efficient generation of large ECG basis sets is essentially impossible when using the `G` option.
- `I` selects the iterative solver based on the inverse iteration method. When `I` is used, the solver relies on the header values `CURRENT_ENERGY`, `EIGVAL_TOLERANCE`, and `INVITPARAMETER` (it targets the eigenvalue close to `CURRENT_ENERGY` $\times$ `INVITPARAMETER`), and `WHICH_EIGENVALUE` is not referenced. The `I` option is much faster than `G`. When it comes to updating the eigenvector/eigenvalue routinely (as is done in `BASIS_ENL`, `OPT_CYCLE`, and `FULL_OPT1`) it may be several orders of magnitude faster. The reason for this is that updating the solution in the inverse iteration approach scales as $\mathcal{O}(K^2)$ , where $K$ is the basis size. Using the `G` option that calls standard LAPACK eigensolver results in $\mathcal{O}(K^3)$ scaling.

This eigenvalue solver type argument is not repeated in detail for each command below.

#### 1. `BASIS_ENL`

Grows the basis by stochastic selection of new basis functions followed by optimization of their nonlinear parameters. At each step a set of new candidate functions is generated a certain number of times (`500` in the example below) based on the existing distribution of nonlinear parameters, and only the set that lowers the energy the most is kept and then optimized.

*Example* :

`BASIS_ENL  G    1    5    1  500  200  0.95  3.0`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `1` | integer | Function number from which the basis enlargement starts (normally the current basis size plus one). If the value exceeds the current basis size plus one, it is automatically reset to the current basis size plus one. |
| `5` | integer | The size to which the basis must be grown. |
| `1` | integer | The number of functions that are randomly selected and added to the basis at each enlargement step. Normally it is best to add only one function at a time, but in some special cases one can consider adding more than one function. |
| `500` | integer | The number of random trials for the stochastic selection at each step. |
| `200` | integer | Maximum number of energy evaluations in the optimization of the nonlinear parameters of the best new candidate(s) that follows the stochastic selection. It should not be too small (e.g. 3-5), as little or no progress will be made by the minimizer, nor too large (typically a few hundred is enough), so that time is not wasted when the minimization gets stuck. The proper value depends on how many functions are added at once and on the number of particles in the system. |
| `0.95` | real | Pair overlap threshold. A new basis function is not added if its overlap (absolute value of $S_{ij}$ ) with any other basis function exceeds the threshold. This prevents building nearly linearly dependent basis sets. Setting this to zero or a negative value disables the check. |
| `3.0` | real | Linear coefficient threshold (for coefficients in front of normalized basis functions). A new basis function is not added if the resulting linear coefficient of any function exceeds the threshold by magnitude (because ECG basis functions form a non-orthogonal basis, the linear coefficients can be greater than `1.0`). This prevents near linear dependencies that manifest as several functions having huge coefficients of opposite sign. Setting this to zero or a negative value disables the check. |

#### 2. `OPT_CYCLE`

Performs a cyclic optimization of the current basis, optimizing one or several functions at a time and shifting the optimized window along the basis. This sweep is repeated a requested number of times.

*Example* :

`OPT_CYCLE  G    5    1    5    1    1   10  130  0.95  3.0   10`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `5` | integer | Current basis size (must match the actual basis size). |
| `1` | integer | Function number from which the optimization cycle begins. Normally it should start from function 1. In some specific situations, however, one might want to skip optimizing the first few functions. |
| `5` | integer | Function number at which the optimization cycle ends. To optimize the entire basis, set the begin function to 1 and the end function equal to the current basis size. |
| `1` | integer | Number of functions whose parameters are optimized simultaneously at each step. Normally set to 1, which is usually the most efficient strategy. For some special cases (e.g. a "molecular" system) optimizing the parameters of more than one function at a time may be helpful. |
| `1` | integer | Number of functions by which the optimized window is shifted at each step. Normally set to 1 (or to the number of functions optimized simultaneously). |
| `10` | integer | Number of optimization cycles (sweeps) to perform. For small basis sizes, when cycles are cheap, one may cycle multiple times. For larger basis sets it is usual to do just one cycle each time after enlarging the basis by several functions. |
| `130` | integer | Maximum number of energy evaluations at each step. See the analogous parameter of `BASIS_ENL` for guidance on choosing it. |
| `0.95` | real | Pair overlap threshold for the penalty function. If an overlap $S_{ij}$ exceeds the threshold, a smooth quadratic penalty is added to the objective function. Setting this to a value equal to or greater than `1.0` or smaller than `0.0` disables the penalty. |
| `3.0` | real | Linear coefficient threshold (for coefficients in front of normalized basis functions). The optimized parameters are not stored if any resulting linear coefficient exceeds the threshold by magnitude. Setting this to zero or a negative value disables the check. |
| `10` | integer | Saving frequency: the updated basis is written to the input/output file after each block of this many cycle steps. For large systems, where calculations are more expensive, `1` is a good choice. For small systems (e.g. 3 particles) and small basis sizes a small value such as `1` may result in excessive disk writing; a value of `5` or `10` is generally a good choice then. |

#### 3. `FULL_OPT1`

Performs a full (i.e. simultaneous) optimization of the nonlinear parameters of all basis functions, or of a selected contiguous subset. Note that this command may require a lot of time and a large amount of memory to store the Hessian matrix (approximately $n_v(n_v+1)/2$ real values, where $n_v=K\times N_p$, $K$ is the basis size, and $N_p$ is the number of nonlinear parameters per function). It is therefore not usually efficient for routine optimization of very large bases (e.g. $> 1000$ functions).

*Example* :

`FULL_OPT1  G    5    1    5 999999  0.95  3.0   60   60  none`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `5` | integer | Current basis size (must match the actual basis size). |
| `1` | integer | First function in the subset of functions whose parameters are to be optimized. |
| `5` | integer | Last function in the subset of functions whose parameters are to be optimized. To optimize the entire basis, set this equal to the current basis size. |
| `999999` | integer | Maximum number of energy evaluations. Normally set to some very large number so that the optimization continues until convergence; a smaller value imposes a hard limit after which the optimization stops. Upon exit the input/output file is updated with the best Gaussian parameters found. |
| `0.95` | real | Pair overlap threshold for the penalty function. If an overlap $S_{ij}$ exceeds the threshold, a smooth quadratic penalty is added to the objective function. Setting this to a value equal to or greater than `1.0` or smaller than `0.0` disables the penalty. |
| `3.0` | real | Maximum overlap penalty, which defines the magnitude of the overlap penalty. With this value equal to 1.0, an overlap of $S_{ij}=1.0$ contributes a penalty of 1.0. |
| `60` | integer | Time interval, in seconds, between successive saves of the current best basis in the input/output file (e.g. 60 means every minute). |
| `60` | integer | Time interval, in seconds, between successive saves of the Hessian matrix (e.g. 60 means every minute). Note that the Hessian may require a large amount of storage (see the note above). |
| `none` | string | Name of the file where the Hessian matrix is stored. Any name is acceptable except `none` (also `None`, `NONE`), which disables Hessian storage altogether. Storing the Hessian allows a more efficient restart after a break or failure, since the Hessian approximation does not need to be rebuilt from multiple gradient evaluations; for very small bases this is unnecessary, while for very large bases it may require prohibitive storage. |

#### 4. `EXPC_VALS`

Computes expectation values for the current basis.

*Example* :

`EXPC_VALS  I    5`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `I` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `5` | integer | Current basis size (must match the actual basis size). |

#### 5. `DENSITIES`

Computes densities of particles in the center-of-mass frame as well as pair correlation functions. Currently this is implemented only for the case of `RG_0S`, `RG_1P`, and `CG_0S` basis types. This command is an extension of `EXPC_VALS`: it does everything `EXPC_VALS` does and, in addition, evaluates densities and correlation functions on user-supplied grids. The grid files contain one grid point per line with no blank lines; for $L=0$ each point is a single radius $r\ge 0$, while for $L=1$ each point is two cylindrical coordinates (distance to the $z$-axis $\rho\ge 0$ and the $z$-coordinate). The output files begin with a `#` header line and reproduce the grid columns followed by the computed quantities; only inequivalent functions are written, e.g. for an $S$-state of beryllium that has four identical electrons, $g_1$ and $g_{12}$ for correlation functions, $\rho_1$ and $\rho_2$ for densities are computed.

*Example* :

`DENSITIES  I    5  cf_grid.dat  cf.dat  dens_grid.dat  dens.dat`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `I` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `5` | integer | Current basis size (must match the actual basis size). |
| `cf_grid.dat` | string | Input file with the grid of points at which the pair correlation functions are evaluated. |
| `cf.dat` | string | Output file with the computed pair correlation functions $g_i$ and $g_{ij}$. |
| `dens_grid.dat` | string | Input file with the grid of points at which the particle densities are evaluated (same format as `cf_grid.dat`; the same file may be used for both). |
| `dens.dat` | string | Output file with the computed densities $\rho_i$ in the center-of-mass frame. |

#### 6. `MOMT_DENS`

Computes momentum densities of particles in the center-of-mass frame as well as momentum pair correlation functions. Currently this is implemented only for the case of `RG_0S` basis type. Like `DENSITIES`, it is an extension of `EXPC_VALS`. The grid file format and output file conventions are the same as for `DENSITIES`, but the computed quantities are evaluated in momentum space.

*Example* :

`MOMT_DENS  I    5  mom_cf_grid.dat  mom_cf.dat  mom_dens_grid.dat  mom_dens.dat`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `I` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `5` | integer | Current basis size (must match the actual basis size). |
| `mom_cf_grid.dat` | string | Input file with the grid of points at which the momentum pair correlation functions are evaluated. |
| `mom_cf.dat` | string | Output file with the computed momentum pair correlation functions. |
| `mom_dens_grid.dat` | string | Input file with the grid of points at which the momentum densities are evaluated (same format as `cf_grid.dat`; the same file may be used for both). |
| `mom_dens.dat` | string | Output file with the computed momentum densities in the center-of-mass frame. |

#### 7. `SAVE_FILE`

Saves the current basis into a separate file. Unlike most commands, `SAVE_FILE` has no eigenvalue solver type argument because it just saves a file and does not invoke any eigenvalue solver.

*Example* :

`SAVE_FILE    5  inout_He_1Se-01-00005.txt`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `5` | integer | Current basis size (must match the actual basis size). |
| `inout_He_1Se-01-00005.txt` | string | Name of the file where the basis is stored, written in the current working directory. It is best to use a unique, self-explanatory name (indicating the system, state/term symbol, eigenvalue number, and basis size) rather than `inout.txt` or the name of an existing file, which would be overwritten. |

#### 8. `SAVE_HSWF`

Saves the Hamiltonian and overlap matrices, the eigenvector of linear coefficients, and the full wave function (both linear and nonlinear parameters of the basis functions). In the matrix output files each element is written on a separate line preceded by its two integer indices; the eigenvector is normalized so that $v'Sv=1$. Any of the four output file names may be set to `none` (also `None`, `NONE`) to skip saving that particular quantity. For example, the second example below saves only the wave function.

*Examples* :

`SAVE_HSWF  I    5  H.txt  S.txt  eigvec.txt  wf.txt`

`SAVE_HSWF  I    5  none  none  none  wavefunction.txt`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `I` | character | Eigenvalue solver type, `G` or `I` (see the note above). |
| `5` | integer | Current basis size (must match the actual basis size). |
| `H.txt` | string | Name of the file where the Hamiltonian matrix is stored, or `none` to skip it. |
| `S.txt` | string | Name of the file where the overlap matrix is stored, or `none` to skip it. |
| `eigvec.txt` | string | Name of the file where the eigenvector of linear coefficients is stored (normalized so that $v'Sv=1$), or `none` to skip it. |
| `wf.txt` | string | Name of the file where the full wave function is stored, including both the linear coefficients and the nonlinear parameters of all basis functions, or `none` to skip it. |

#### 9. `ELIM_LCFN`

Eliminates basis functions whose contribution to the energy is small, i.e. those whose linear coefficient (in front of the normalized function) has an absolute value smaller than a given threshold. The reduced basis is written to the specified file and the program then terminates. This command only works with the `G` eigenvalue solver.

*Example* :

`ELIM_LCFN  G    5  1.0E-4  inout_reduced.txt`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type. Must be `G`; the `I` solver is not supported for this command. |
| `5` | integer | Current basis size (must match the actual basis size). |
| `1.0E-4` | real | Linear coefficient threshold. Functions whose linear coefficient (in front of the normalized function) is smaller than this value by magnitude are eliminated. |
| `inout_reduced.txt` | string | Name of the file where the reduced basis is stored. |

#### 10. `ELIM_LND1`

Eliminates linearly dependent functions. It checks pair linear dependency only, removing each function whose overlap (absolute value) with any earlier (smaller-numbered) function exceeds the threshold. The reduced basis is written to the specified file and the program then terminates. This command only works with the `G` eigenvalue solver.

*Example* :

`ELIM_LND1  G    5  0.99  inout_reduced.txt`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type. Must be `G`; the `I` solver is not supported for this command. |
| `5` | integer | Current basis size (must match the actual basis size). |
| `0.99` | real | Pair linear dependency threshold. A function is removed if its overlap (absolute value, computed with normalized functions) with any earlier function exceeds this value. |
| `inout_reduced.txt` | string | Name of the file where the reduced basis is stored. |

#### 11. `SEPR_LND1`

Does the same linear dependency detection as `ELIM_LND1`, but instead of discarding the offending functions it randomly perturbs their nonlinear parameters to separate them. Each affected parameter $a_{\mathrm{old}}$ is replaced by a random value drawn from the interval $[a_{\mathrm{old}}(1-s),\,a_{\mathrm{old}}(1+s)]$, where $s$ is the separation parameter. The result is written to the specified file and the program then terminates. This command only works with the `G` eigenvalue solver.

*Example* :

`SEPR_LND1  G    5  0.99  0.1  inout_separated.txt`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type. Must be `G`; the `I` solver is not supported for this command. |
| `5` | integer | Current basis size (must match the actual basis size). |
| `0.99` | real | Linear dependency threshold (see `ELIM_LND1`). Functions whose overlap exceeds this value are separated rather than removed. |
| `0.1` | real | Separation parameter $s$ controlling the random shift of the nonlinear parameters of the affected functions. |
| `inout_separated.txt` | string | Name of the file where the resulting basis is stored. |

#### 12. `SEPR_FLCF`

Randomly perturbs the nonlinear parameters of basis functions whose linear coefficient (in front of the normalized function) exceeds the threshold by magnitude, in order to separate near linear dependencies that manifest as large coefficients. Each affected parameter $a_{\mathrm{old}}$ is replaced by a random value from the interval $[a_{\mathrm{old}}(1-s),\,a_{\mathrm{old}}(1+s)]$. The result is written to the specified file and the program then terminates. This command only works with the `G` eigenvalue solver.

*Example* :

`SEPR_FLCF  G    5  3.0  0.1  inout_separated.txt`

*Arguments and their description* :

| Argument | Type | Description |
| :--- | :---: | :--- |
| `G` | character | Eigenvalue solver type. Must be `G`; the `I` solver is not supported for this command. |
| `5` | integer | Current basis size (must match the actual basis size). |
| `3.0` | real | Linear coefficient threshold. Functions whose linear coefficient (in front of the normalized function) exceeds this value by magnitude are separated. |
| `0.1` | real | Separation parameter $s$ controlling the random shift of the nonlinear parameters of the affected functions. |
| `inout_separated.txt` | string | Name of the file where the resulting basis is stored. |

### History

The history section is generated and updated by the program as the basis is being built; the user does not write it when starting a calculation from scratch. It is surrounded by two ruler lines (a row of `=` characters) and contains one line for every basis size from 1 up to the current basis size. Each line records the lowest variational energy obtained at the corresponding basis size, together with three integer counters that allow the program to resume a calculation that was interrupted without repeating work that has already been done (see also the section on restarting calculations below). The last line of the history (the one for the current basis size) is the same line that is reproduced in the header excerpt.

A typical history section looks as follows (this example corresponds to the ten-function He ground-state basis from the example at the top of this page):

```
 ==============================
      1 -0.1668865998401775E+01      0      0      0
      2 -0.2642629552717170E+01      0      0      0
      3 -0.2798648895724608E+01      0      0      0
      4 -0.2844201302182638E+01      0      0      0
      5 -0.2889258918771523E+01     10      0    121
      6 -0.2893514456170832E+01      0      0      0
      7 -0.2895565511375141E+01      0      0      0
      8 -0.2896518341033143E+01      0      0      0
      9 -0.2897179931311348E+01      0      0      0
     10 -0.2900744334456869E+01     10      0    123
 ==============================
```

The columns, from left to right, are:

| Column | Type | Description |
| :--- | :---: | :--- |
| 1 | integer | Basis size, i.e. the number of basis functions to which this line refers. The lines are listed in order of increasing basis size, from 1 to the current basis size. |
| 2 | real | The lowest variational energy (for the targeted eigenvalue) that was obtained when the basis had this number of functions. |
| 3 | integer | The number of cyclic optimization cycles (`OPT_CYCLE`) that have already been completed at this basis size. It is used to resume a cyclic optimization at the correct cycle if the calculation is restarted. |
| 4 | integer | The function number at which the last completed optimization step began at this basis size. It is used to resume a cyclic optimization at the correct function within a cycle if the calculation is restarted. |
| 5 | integer | The number of energy evaluations already spent on full optimization (`FULL_OPT1`) at this basis size. It is used to enforce the maximum-number-of-energy-evaluations limit across restarts. |

All three integer counters are reset to zero whenever a new function is added at a given basis size (e.g. by `BASIS_ENL`), and a value of zero therefore indicates that no cyclic or full optimization has yet been performed at that basis size.

### Basis functions

The basis functions section is the last section of the input file. Like the history, it is generated and updated by the program and is not written by the user when a calculation is started from scratch. It contains the nonlinear variational parameters of the basis functions, one function per line, listed in order of increasing function index. Note that the linear expansion coefficients are not stored here: they are not independent variational parameters but are obtained by solving the generalized eigenvalue problem, and are therefore recomputed whenever the basis is read.

For the `RG_0S` basis type, each line consists of the integer function index followed by the $n_{pt}=n(n+1)/2$ nonlinear parameters of that function, where $n$ is the number of pseudoparticles (one less than the number of particles). These parameters are the independent elements of the lower-triangular matrix that defines the exponent of the explicitly correlated Gaussian. A typical section for the ten-function He ground-state basis from the example at the top of this page (where $n=2$, so each function has three nonlinear parameters) looks as follows:

```
 ==============================
      1  0.8225526505330379E+00  0.3286331967568486E-01  0.4734092048850279E+00
      2  0.4657831516551231E+00 -0.4035888660117890E-01  0.1841499686869463E+01
      3  0.5792314542342561E+00 -0.2455636629377815E-01  0.8488492583438731E+01
      4  0.1316320598387222E+01  0.1264307213755476E+00  0.7672315752278411E+00
      5  0.1230152783754290E+01  0.1911586281826508E+00  0.2384612056831933E+01
      6 -0.2233523534916765E+01  0.1516523349832122E+01 -0.1098568271679442E+01
      7  0.9104134991273443E+00  0.2902311989845964E-01  0.1975856897728040E+02
      8  0.1405523968799062E+01  0.3635708756466438E+00  0.6425523571898003E+00
      9 -0.6106458148603032E+00 -0.6976991906481144E-02 -0.3667310281682786E+01
     10  0.1577666796470805E+01  0.2898602080507635E+00  0.5856740074044248E+01
```

The other basis types insert one or two integer indices between the function index and the nonlinear parameters: the `RG_1P` basis type has one additional integer (the $z$-index), while the `RG_2D` and `RG_2P` basis types have two additional integers. The number of nonlinear parameters per function, $n_{pt}$, and their precise meaning depend on the basis type.

The complex Gaussian basis `CG_0S` has twice as many Gaussian nonlinear parameters (for the same number of pseudoparticles $n$ ) compared to the real bases `RG_0S`, `RG_1P`, `RG_2D` and `RG_2P`.

## Creating input files

When a large basis needs to be generated it requires multiple actions (e.g. basis enlargements `BASIS_ENL` and cyclic optimizations `OPT_CYCLE`) that need to be repeated routinely. It may be tedious to manually write the sequence of corresponding commands and their multiple arguments. For this reason it is convenient to create a simple utility that creates a desired list of commands, or, even better, prepares a ready-to-use input file `inout.txt` to initiate calculations from scratch.

Directory `ecgpack/utilities/input_file_manipulation` contains a Python script `ecg_input_file_generator.py` that can be invoked to generate an input file for a desired system with a long list of repeated actions. One can modify this script as needed.

## Restarting calculations that generate or optimize a basis

As the program saves/overwrites the input file `inout.txt` regularly when the basis is extended or optimized, restarting calculations from the point where they were terminated is trivial and automatic. The only exception may take place when one invokes the full optimization (command `FULL_OPT1`). A proper restart of the full optimization requires regular saving of the Hessian into a file. The user can do this by replacing `none` (which means no Hessian file saving) with the name of the file, e.g.

`FULL_OPT1  G    5    1    5 999999  0.95  3.0   60   60  hessian.dat`

It should be kept in mind, however, that saving a Hessian in the case of large basis sets may take an enormous amount of disk space. So this option should be used with caution.

If the Hessian was not saved, the full optimization will still resume, but it may not proceed as efficiently as it was proceeding before termination because the information about the Hessian approximation that was built previously was lost.
