# libtsm

A Python toolbox for Time-Scale Modification (TSM) and Pitch-Shifting.

Details and example application:

https://www.audiolabs-erlangen.de/resources/MIR/2021-DAFX-AdaptivePitchShifting

`libtsm` is based on a re-implementation of the Matlab TSM Toolbox by Jonathan Driedger and Meinard Müller:

https://www.audiolabs-erlangen.de/resources/MIR/TSMtoolbox/

If you use the libtsm in your research, please consider the following references.

## References

Sebastian Rosenzweig, Simon Schwär, Jonathan Driedger, and Meinard Müller:
Adaptive Pitch-Shifting with Applications to Intonation Adjustment in A Cappella Recordings
Proceedings of the International Conference on Digital Audio Effects (DAFx), 2021.

Jonathan Driedger and Meinard Müller:
[TSM Toolbox: MATLAB Implementations of Time-Scale Modification Algorithms.](https://www.audiolabs-erlangen.de/fau/professor/mueller/publications/2014_DriedgerMueller_TSM-Toolbox_DAFX.pdf)
In Proceedings of the International Conference on Digital Audio Effects (DAFx): 249–256, 2014.

Jonathan Driedger and Meinard Müller:
[A Review on Time-Scale Modification of Music Signals.](https://www.mdpi.com/2076-3417/6/2/57)
Applied Sciences, 6(2): 57–82, 2016.

Jonathan Driedger, Meinard Müller, and Sebastian Ewert:
[Improving Time-Scale Modification of Music Signals using Harmonic-Percussive Separation.](https://ieeexplore.ieee.org/abstract/document/6678724)
IEEE Signal Processing Letters, 21(1): 105–109, 2014.

## Installation

With Python >= 3.6, you can install libtsm using the Python package manager pip:

```
pip install libtsm
```

## Documentation

The API documentation of `libtsm` is hosted here:

https://meinardmueller.github.io/libtsm

## Contributing

We are happy for suggestions and contributions. However, to facilitate the synchronization, we would be grateful for either directly contacting us via email (meinard.mueller@audiolabs-erlangen.de) or for creating [an issue](https://github.com/meinardmueller/libtsm/issues) in our GitHub repository. Please do not submit a pull request without prior consultation with us.

If you want to report an issue with libtsm or seek support, please use the same communication channels (email or GitHub issue).

## Tests

Central to our tests is the comparison of `libtsm` with the MATLAB TSM Toolbox.
To this end, please execute `tests/test_matlab.m` in MATLAB to create the MATLAB output.
Then, you can use [pytest](https://pytest.org) for executing our Python test scripts. `pytest` is available when installing libtsm with the extra requirements for testing.

```
pip install 'libtsm[tests]'
pytest
```

## Acknowledgements

This project is supported by the German Research Foundation (DFG MU 2686/12-1, MU 2686/13-1).
The International Audio Laboratories Erlangen are a joint institution of the Friedrich-Alexander Universität Erlangen-Nürnberg (FAU) and Fraunhofer Institute for Integrated Circuits IIS. We thank Edgar Suarez, El Mehdi Lemnaouar and Miguel Gonzales for implementation support.
