#####################
# Configuration.pm

package Configuration;
#-# use lib "/home/cgrb/cgrblib-dev/perl5/Diurnal";
use FindBin;
use lib "$FindBin::Bin/../Library";
use lib '/deepspace/mockler/priesth/Scripts';
use Config::Tiny;

# Designed to read a .ini formatted configuration file.
# Example of .ini:

# [Section]
# name=value
# name=value
# name=value

# [Section2]
# name=value
# name=value
# etc..

# This module currently is basically just a wrapper for Config::Tiny.  The reason
# I built this wrapper was so that in the future if another Config module needed to be used,
# we can simply replace references to Config::Tiny here with the new module, and change the 
# syntax of how we obtain the values from the .ini formated config file.

sub new
{
	# allow for inheritance
	my $class=shift;
	my $this            		= {};
	my $file=shift;
	$this->{config}			= Config::Tiny->new();
	$this->{config}			= Config::Tiny->read($file);
	
	# allow for inheritance, if needed
	bless($this,$class);
	return $this;
} # end sub new

sub getAll
{
	my $this=shift;
	my $key=shift;
	return keys %{$this->{config}{$key}}
}

# $config->get('SECTION','variable');
sub get
{
	my ($this,$section,$key)	 	= @_;
	if(exists $this->{config}->{$section}->{$key} ) { return $this->{config}->{$section}->{$key}; }
	return -1;
}



###########################################################################
# 			Perldoc begin
###########################################################################

=head1 Name

Config.pm - Used to extract values from a .ini formatted configuration file.

=head1 Synopsis

# configuration file that is .ini formatted like this:

# [MYSQL]
# host=localhost

# perl code to access .ini formatting
require Config;

my $config	= new Config(file=>'config.conf');
my $host	= $config->get("MYSQL","host"); # get(section,variable)

=head1 Description

This is a simple module to deal with simple configuration files that are .ini formatted.
Currently, this module is a wrapper to Config::Tiny (cpan).  I wrapped the module so that in
the future if there is ever a need to change configuration modules (for something more complicated, 
or just different), the code will only have to be changed in one place instead of everywhere 
information from the configuration file is used.

=head1 Author

Adam Gustafson

=cut







1;

