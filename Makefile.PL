use ExtUtils::MakeMaker;
WriteMakefile(
    NAME            => 'Bio::NGS::HCL::DeCovA',
    ABSTRACT        => "To be filled...",
    #VERSION_FROM    => 'lib/Your/Module.pm'
    VERSION_FROM    => 'bin/DeCovA',
    AUTHOR          => [
	"Thomas Simonet <Thomas.Simonet@chu-lyon.fr>",
	"Pierre-Antoine Rollat-Farnier <Pierre-Antoine.Rollat-Farnier@chu-lyon.fr>",
    ],
    EXE_FILES       => [ "bin/DeCovA" ],
);

