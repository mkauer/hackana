#! /bin/sh
#
# olbdd	       Start/Stop the OLBD daemon
#
# chkconfig: 345 20 80
# description:	The olbd daemon runs the XRD control network.
#
# processname: olbd
# pidfile: /var/run/olbd.pid
# config:

OLBD=$(ROOTSYS)/bin/olbd
OLBDLIBS=@libdir@

# Source function library.
. /etc/init.d/functions

# Get config.
. /etc/sysconfig/network

# Get olbd config
[ -f /etc/sysconfig/olbd ] && . /etc/sysconfig/olbd

# Configure the environment
[ ! -z "$OLBDENVCONFIG" ] && [ -f "$OLBDENVCONFIG" ] && . $OLBDENVCONFIG

# Check that networking is up.
if [ ${NETWORKING} = "no" ]
then
	exit 0
fi

[ -x $OLBD ] || exit 0

RETVAL=0
prog="olbd"

start() {
        echo -n $"Starting $prog: "
        # Options are specified in /etc/sysconfig/olbd .
        # See $ROOTSYS/etc/daemons/olbd.sysconfig for an example.
        # $OLBDUSER *must* be the name of an existing non-privileged user.
        export LD_LIBRARY_PATH=$OLBDLIBS:$LD_LIBRARY_PATH
        daemon $OLBD -b -l $OLBDLOG -R $OLBDUSER -c $OLBDCF $OLBDDEBUG
        RETVAL=$?
        echo
        [ $RETVAL -eq 0 ] && touch /var/lock/subsys/olbd
        return $RETVAL
}

stop() {
	[ ! -f /var/lock/subsys/olbd ] && return 0 || true 
        echo -n $"Stopping $prog: "
        killproc olbd
        RETVAL=$?
        echo
        [ $RETVAL -eq 0 ] && rm -f /var/lock/subsys/olbd
	return $RETVAL
}

# See how we were called.
case "$1" in
  start)
	start
	;;
  stop)
	stop
	;;
  status)
	status olbd
	RETVAL=$?
	;;
  restart|reload)
	stop
	start
	;;
  condrestart)
	if [ -f /var/lock/subsys/olbd ]; then
            stop
            start
        fi
	;;
  *)
	echo  $"Usage: $0 {start|stop|status|restart|reload|condrestart}"
	exit 1
esac

exit $RETVAL
