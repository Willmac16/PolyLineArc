import numpy as np, math
import matplotlib.pyplot as plt

TANGENCY_THRESHOLD = 1e-13
COINCIDENT_THRESHOLD = 1e-16

def offset_vector(heading, offset):
    clockwise = np.array(((0, 1),(-1, 0)))
    vec = np.matmul(clockwise, heading)
    vec *= offset / np.linalg.norm(vec)
    return vec

class PolyLineArc:
    def __init__(self, next):
        if next != None:
            if np.linalg.norm(next.start - self.end) < COINCIDENT_THRESHOLD:
                self.next = next
            else:
                raise Exception("PolyLineArc: next segment must start at end of current segment")
        return

    # Offset is a signed float ninety degrees to the right of heading
    def offset(self, offset: float):

        end_offset_vector = offset_vector(self.end_heading(), offset)

        if self.next is None:
            return self.offset_self(None, offset)

        # If two segments aren't tangent figure out if we will procede
        if np.linalg.norm(self.end_heading() - self.next.start_heading()) >= TANGENCY_THRESHOLD:
            # Check if the two segments need a fillet patch or an intersection blend
            heading_delta = self.next.start_heading() - self.end_heading()

            # print(heading_delta)
            if np.dot(heading_delta, end_offset_vector) >= TANGENCY_THRESHOLD:
                # intersection blend
                # nope I'm just gonna throw an error
                raise Exception("Illegal convex a-tangent intersection")
                # Convex intersections need to be tangent cux they aren't machinable
                # and I'm not about to do arc arc intersection code
            else:
                o_next = self.next.offset(offset)

                # fillet patch
                patch_start = self.end + end_offset_vector
                patch_end = self.next.start + offset_vector(self.next.start_heading(), offset)
                patch = PolyArc(patch_start, patch_end, self.end, o_next)

                return self.offset_self(patch, offset)
        else:
            o_next = self.next.offset(offset)
            return self.offset_self(o_next, offset)

    def length(self):
        return

    def pos(self, t):
        return

    def pos_prime(self, t):
        return

    def offset_self(self, next, offset):
        return

    def start_heading(self):
        return

    def end_heading(self):
        return

    def plot(self):
        xs, ys = self.plot_vals()

        plt.plot(xs, ys)
        plt.title("PolyLineArc")

    def plot_vals(self):
        xs = []
        ys = []

        for t in np.arange(0, 1.01, 0.01):
            pos = self.pos(t)
            xs.append(pos[0])
            ys.append(pos[1])

        if self.next is not None:
            n_xs, n_ys = self.next.plot_vals()

            xs.extend(n_xs)
            ys.extend(n_ys)

        return (xs, ys)

    def plot_prime(self):
        xs, ys = self.plot_vals_prime()

        plt.plot(xs, ys)
        plt.title("PolyLineArc Slope")

    def plot_vals_prime(self):
        xs = []
        ys = []

        for t in np.arange(0, 1.01, 0.01):
            head = self.pos_prime(t)
            pos = self.pos(t)
            xs.append(pos[0])
            ys.append(head[1]/head[0])

        if self.next is not None:
            n_xs, n_ys = self.next.plot_vals_prime()

            xs.extend(n_xs)
            ys.extend(n_ys)

        return (xs, ys)

    def num_segments(self):
        if self.next is None:
            return 1
        else:
            return 1 + self.next.num_segments()

class PolyLine(PolyLineArc):
    def __init__(self, start: np.ndarray, end: np.ndarray, next: PolyLineArc):
        self.start = start
        self.end = end

        self.next = None

        super().__init__(next)

    def length(self):
        delta = self.end - self.start

        return np.linalg.norm(delta)

    def pos(self, t):
        return self.start + (self.end - self.start) * t

    def pos_prime(self, t):
        return self.end - self.start

    def offset_self(self, next, offset):
        radial = offset_vector(self.heading(), offset)

        off_line = PolyLine(self.start + radial, self.end + radial, next)
        return off_line

    def heading(self):
        connect_vec = self.end - self.start
        connect_vec /= np.linalg.norm(connect_vec)
        return connect_vec

    def start_heading(self):
        return self.heading()

    def end_heading(self):
        return self.heading()

    def fillet(self, fillet_radius):
        if self.next:
            self.next.fillet(fillet_radius)

            if self.next.heading != self.heading:
                intersection = self.end

                u = self.start - intersection
                v = self.next.end - intersection

                u /= np.linalg.norm(u)
                v /= np.linalg.norm(v)

                angle = math.acos(np.dot(u, v))
                trim_length = fillet_radius * math.tan((math.pi - angle) / 2)

                new_end = intersection + trim_length * u
                new_next_start = intersection + trim_length * v

                self.end = new_end
                self.next.start = new_next_start

                fillet_arc = PolyArc(new_end, new_next_start, np.zeros(1), self.next, radius=fillet_radius)

                self.next = fillet_arc

                # print("line to fillet", np.equal(self.end, self.next.start).all())
                # print("fillet to line", np.equal(self.next.end, self.next.next.start).all())




class PolyArc(PolyLineArc):
    # If radius is defined then ignore center and solve for center given start and end
    def __init__(self, start: np.ndarray, end: np.ndarray, center: np.ndarray, next: PolyLineArc, radius: float = None):
        self.start = start
        self.end = end

        if radius == None:
            self.center = center
        else:
            # Code to Find the center of the arc given the start, end, and diameter
            # I am just stealing this from my prior work with gcode G2 R commands
            directConnect = end - start

            q = np.array((directConnect[1]*-1, directConnect[0]))
            q /= np.linalg.norm(q)
            q *= -1 * math.sqrt(radius**2 - (np.linalg.norm(directConnect)/2)**2)

            self.center = start + directConnect/2 + q

        self.next = None

        super().__init__(next)

    def angle(self):
        r_start = np.append(self.start_heading(), (0))
        r_end = np.append(self.end_heading(), (0))

        return math.asin(np.cross(r_start, r_end)[2])

    def length(self):
        return self.radius() * abs(self.angle())

    def pos(self, t):
        start_radius = self.start - self.center

        angle = self.angle()

        working_angle = angle * t

        rot_matrix = np.array(((np.cos(working_angle), -np.sin(working_angle)), (np.sin(working_angle), np.cos(working_angle))))

        return self.center + np.matmul(rot_matrix, start_radius)

    def pos_prime(self, t):
        working_angle = self.angle() * t

        rot_matrix = np.array(((np.cos(working_angle), -np.sin(working_angle)), (np.sin(working_angle), np.cos(working_angle))))

        return np.matmul(rot_matrix, self.start_heading())

    # positive for anti-clock, negative for clock
    def turn_handedness(self):
        r_start = self.start - self.center
        r_end = self.end - self.center

        return (np.cross(r_start, r_end) > 0)*2 - 1

    def offset_self(self, next, offset):
        radial_start = offset_vector(self.start_heading(), offset)
        radial_end = offset_vector(self.end_heading(), offset)

        off_arc = PolyArc(self.start + radial_start, self.end + radial_end, self.center, next)
        return off_arc

    def start_heading(self):
        r_start = self.start - self.center
        r_start /= np.linalg.norm(r_start)

        r_start *= self.turn_handedness()

        rotate = np.array(((0, -1),(1, 0)))

        return np.matmul(rotate, r_start)

    def end_heading(self):
        r_end = self.end - self.center
        r_end /= np.linalg.norm(r_end)

        r_end *= self.turn_handedness()

        rotate = np.array(((0, -1),(1, 0)))

        return np.matmul(rotate, r_end)

    def radius(self):
        return np.linalg.norm(self.center - self.start)



def test():

    arc = PolyArc(np.array((-0.018604, -0.015926)), np.array((0.005715, -0.005080)), np.array((0.050201, -0.137508)), None)
    arc_2 = PolyArc(np.array((-0.018604, -0.015926)), np.array((0.005715, -0.005080)), np.array((0.0)), None, radius=0.1)

    print(arc.radius())
    print(arc_2.radius())

    arc.plot()
    arc_2.plot()
    plt.show()



if __name__ == "__main__":
    test()
