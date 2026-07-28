# Writing the wave function

## Description

The sample `inout.txt` file in this directory defines a basis of 100 ECG functions for the lowest $^2S^e$ state of the lithium atom (assuming infinite nuclear mass). Running it generates the wave-function expansion and writes it to a separate file whose name is specified in the input (`wavefunction.txt`).

On a single CPU core, the calculation takes approximately 1–2 seconds in double precision (gfortran compiler, AMD Ryzen 7 7800X3D).

## The `SAVE_HSWF` instruction

To produce the wave function, the following instruction is added to `inout.txt`:

```text
SAVE_HSWF I 100 none none none wavefunction.txt
```

The relevant arguments are:

| Argument | Description |
|---|---|
| `SAVE_HSWF` | Keyword for writing the wave function |
| `I` | Eigensolver type |
| `100` | Number of ECG basis functions |
| `wavefunction.txt` | Name of the saved wave-function file |

The basis-set size given in the `SAVE_HSWF` instruction must match the size of the wave function being saved.