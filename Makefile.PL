use ExtUtils::MakeMaker;
WriteMakefile(
    NAME            => 'DeCovA',
    ABSTRACT        => "depth and coverage analysis and CNV detection for NGS data",
    #VERSION_FROM    => 'lib/Your/Module.pm'
    VERSION_FROM    => 'bin/DeCovA',
    AUTHOR          => [
	"Thomas Simonet <Thomas.Simonet@chu-lyon.fr>",
	"Pierre-Antoine Rollat-Farnier <Pierre-Antoine.Rollat-Farnier@chu-lyon.fr>",
    ],
    EXE_FILES       => [ "bin/DeCovA" ],
);

