# Finding-Connections

The communication network at a given bank is composed of several departmental networks, each of which connects to the remaining networks only through a server. This
network topology was chosen for security reasons. All critical data sent between
departmental networks are stored on the hard disk of these servers in order to be
analysed by automatic filtering tools. In case of a security breach at a given department,
the contents of the server's hard disk are instantaneously moved to another
departmental server and the local network becomes disconnected from the remaining
networks.
<br><br>
Unfortunately, the current network cables do not support fast copies between servers. In
order to improve this situation, the IT department is analysing several possibilities to
add a faster cable that only connects these servers. One option is to directly connect all
pairs of servers by cable, which is extremely fast, but expensive. The other option is to
use a tree topology, which is cheaper but also slower.
<br>
This is where you come in. You are doing an intership at the IT department and you have
to implement a program that, given a map of the complete communication network at
the Bank, it locates all departmental servers and computes the amount of cable that is
required for both options. In order to be cost effective, you have to consider the least
amount of cable and use the current infrastruture, that is, the new cable will go through
the buildings side-by-side with the old cable, possibly passing near other network
equipment.
<br><br>
Since you do not know what to expect from the data, your program should be prepared
to handle strange cases: there is no server, or there a single server and therefore, there
is no need to buy new cable; some parts of the network may be completely isolated, but
they may contain servers, for which the new cable should also be installed among them -
to cut expenses, no connection with servers from the main communication network is
required.

# Input

The first line of each test case gives the amount of network equipment (n ≤ 1000) at the
communication network. Note that some nodes of this network may be servers, but you
do not know which ones. Then, the next lines provide information about the network
infrastructure. Each line contains a pair of positive integers that correspond to the
internal id of the network equipment that is directly connected with the existing cable,
followed by the length of this cable, as a positive integer value.
Each test case terminates with number 0. Then, other test cases may follow.

# Output

For each test case, print the number of servers, the total amount of cable for the case of
a fully connected network and the total amount of cable for a tree topology, as required
by the IT department. If there is no server, report "no server".
