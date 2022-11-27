project(
	'neolithic_cpp',
	'cpp',
	version: '0.1.0',
	license: 'GPL',
	meson_version: '>=0.47',
	default_options: ['cpp_std=c++17', 'warning_level=2'],
)

cairo = dependency('cairo', required: get_option('use_cairo'))
omp = dependency('openmp', required: get_option('use_openmp'))
threads = dependency('threads')

cmdeps = []
n2wdeps = [threads]
if omp.found() and not get_option('use_openmp').disabled()
	cmdeps = [omp]
	n2wdeps = [omp, threads]
endif


create_matrix = executable('cm', ['create_helper_matrix.cpp', 'read_table_cpp.cpp'], dependencies: [omp])
neolithic2w = executable('n2w', ['neolithic2w.cpp', 'read_table_cpp.cpp'], dependencies: [omp, threads])
neolithic2 = executable('n2', ['neolithic2.cpp', 'read_table_cpp.cpp'], dependencies: [omp, threads])
spas = executable('spas', ['spaggr_simple.cpp', 'read_table_cpp.cpp'])
sr = executable('sr', ['sample_res.cpp', 'read_table_cpp.cpp'])
if cairo.found() and not get_option('use_cairo').disabled()
	sp2png = executable('sp2png', ['spout2png.cpp', 'read_table_cpp.cpp', 'cellnet_img.cpp'], dependencies: [cairo])
endif
rm = executable('rm', 'runmultiple.cpp', dependencies: [threads])
