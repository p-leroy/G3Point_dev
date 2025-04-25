%% compare matlab:ellipsoid_ex2im and python:explicit_to_implicit

center = [0.21357835, 5.81009748, 0.16241896];
radii = [0.21137105, 0.07781998, 0.0815587 ];
R = [-0.9263723  -0.36679569  0.0854124;  ...
       -0.30781516  0.86809616  0.38943406; ...
       0.21698891 -0.33446969  0.91708551];

matlab_p = ellipsoid_ex2im(center, radii, R);

% for the stack 0, python, we obtain the following implicit parameters
python_p = [   41.93209985,   144.26738821,   151.64467324,   -94.85891829, ...
          16.70173584,    18.01863136,   530.51530604, -1659.08193363, ...
        -157.51707559,  4774.85246949];

%%
matlab_p' - python_p

%%
[center, radii, quat, R] = ellipsoid_im2ex(python_p);

%%
plot_ellipsoid_im(python_p);
