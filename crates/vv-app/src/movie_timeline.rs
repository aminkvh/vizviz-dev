//! The Movie panel's time ruler and track lanes: clips as bars you drag
//! to move and pull at the ends to resize, and a playhead you scrub.

use egui::{
    Align, Align2, Color32, CornerRadius, FontId, Id, Layout, Pos2, Rect, Response, Sense, Stroke,
    StrokeKind, Ui, UiBuilder, Vec2,
};
use egui_phosphor::regular as icon;

use crate::movie::{Clip, Grab, TrackKind};
use crate::theme::{radius, space, text, Tokens};
use crate::ui::AppUi;
use crate::widgets;

const LABEL_WIDTH: f32 = 168.0;
const ROW_HEIGHT: f32 = 28.0;
const RULER_HEIGHT: f32 = 20.0;
const HANDLE_WIDTH: f32 = 6.0;
const MIN_TICK_PX: f32 = 56.0;
const TICK_STEPS: [f32; 12] = [
    0.1, 0.25, 0.5, 1.0, 2.0, 5.0, 10.0, 15.0, 30.0, 60.0, 120.0, 300.0,
];

/// Where seconds sit on screen: the lane area to the right of the labels.
#[derive(Clone, Copy)]
struct TimeAxis {
    left: f32,
    pixels_per_second: f32,
}

impl TimeAxis {
    fn new(area: Rect, duration: f32) -> Self {
        Self {
            left: area.left(),
            pixels_per_second: area.width() / duration.max(0.1),
        }
    }

    fn x(self, seconds: f32) -> f32 {
        self.left + seconds * self.pixels_per_second
    }

    fn seconds(self, x: f32) -> f32 {
        (x - self.left) / self.pixels_per_second
    }
}

/// The smallest tick spacing, in seconds, that keeps labels apart.
fn tick_step(pixels_per_second: f32) -> f32 {
    TICK_STEPS
        .into_iter()
        .find(|s| s * pixels_per_second >= MIN_TICK_PX)
        .unwrap_or(TICK_STEPS[TICK_STEPS.len() - 1])
}

fn clip_color(kind: TrackKind, t: &Tokens) -> Color32 {
    match kind {
        TrackKind::Rotate => t.primary,
        TrackKind::Zoom => t.success,
        TrackKind::MoveTo => t.warning,
        TrackKind::Play => t.text_muted,
    }
}

impl AppUi<'_> {
    pub(crate) fn movie_timeline(&mut self, ui: &mut Ui) {
        let rows = self.movie.movie.tracks.len();
        let size = Vec2::new(
            ui.available_width(),
            RULER_HEIGHT + rows as f32 * ROW_HEIGHT,
        );
        let (rect, _) = ui.allocate_exact_size(size, Sense::hover());
        let lanes = Rect::from_min_max(
            Pos2::new(rect.left() + LABEL_WIDTH, rect.top()),
            Pos2::new(rect.right() - space::GAP, rect.bottom()),
        );
        let axis = TimeAxis::new(lanes, self.movie.movie.duration);
        self.draw_ruler(ui, lanes, axis);
        for track in 0..rows {
            let top = rect.top() + RULER_HEIGHT + track as f32 * ROW_HEIGHT;
            let row = Rect::from_min_size(
                Pos2::new(rect.left(), top),
                Vec2::new(rect.width(), ROW_HEIGHT),
            );
            self.track_label(ui, track, row);
            self.track_lane(ui, track, row, axis);
        }
        self.draw_playhead(ui, rect, axis);
    }

    fn draw_ruler(&mut self, ui: &mut Ui, lanes: Rect, axis: TimeAxis) {
        let t = Tokens::current(ui.ctx());
        let ruler = Rect::from_min_size(lanes.min, Vec2::new(lanes.width(), RULER_HEIGHT));
        let response = ui.interact(ruler, Id::new("movie-ruler"), Sense::click_and_drag());
        if let (true, Some(pos)) = (
            response.clicked() || response.dragged(),
            response.interact_pointer_pos(),
        ) {
            self.movie.playing = false;
            self.movie.set_time(axis.seconds(pos.x));
        }
        let painter = ui.painter();
        painter.rect_filled(ruler, CornerRadius::same(radius::CONTROL), t.surface_raised);
        let step = tick_step(axis.pixels_per_second);
        let font = FontId::proportional(text::CAPTION);
        let ticks = (self.movie.movie.duration / step).floor() as usize;
        for n in 0..=ticks {
            let (seconds, x) = (n as f32 * step, axis.x(n as f32 * step));
            painter.line_segment(
                [
                    Pos2::new(x, ruler.bottom() - 6.0),
                    Pos2::new(x, ruler.bottom()),
                ],
                Stroke::new(1.0, t.text_muted),
            );
            let label = if step < 1.0 {
                format!("{seconds:.1}s")
            } else {
                format!("{seconds:.0}s")
            };
            painter.text(
                Pos2::new(x + 3.0, ruler.top() + 2.0),
                Align2::LEFT_TOP,
                label,
                font.clone(),
                t.text_muted,
            );
        }
    }

    fn draw_playhead(&self, ui: &Ui, rect: Rect, axis: TimeAxis) {
        let t = Tokens::current(ui.ctx());
        let x = axis.x(self.movie.time);
        ui.painter().line_segment(
            [Pos2::new(x, rect.top()), Pos2::new(x, rect.bottom())],
            Stroke::new(2.0, t.danger),
        );
    }

    /// The track's name, and its add-clip and remove buttons.
    fn track_label(&mut self, ui: &mut Ui, track: usize, row: Rect) {
        let label = Rect::from_min_size(row.min, Vec2::new(LABEL_WIDTH - space::GAP, ROW_HEIGHT));
        let kind = self.movie.movie.tracks[track].kind;
        let (mut add, mut remove) = (false, false);
        ui.scope_builder(
            UiBuilder::new()
                .max_rect(label)
                .layout(Layout::left_to_right(Align::Center)),
            |ui| {
                ui.label(format!("{}  {}", track + 1, kind.label()));
                ui.with_layout(Layout::right_to_left(Align::Center), |ui| {
                    remove = widgets::icon_button(ui, icon::TRASH, false)
                        .on_hover_text("Remove this track")
                        .clicked();
                    add = widgets::icon_button(ui, icon::PLUS, false)
                        .on_hover_text("Add a clip at the playhead")
                        .clicked();
                });
            },
        );
        if add {
            self.add_clip_at_playhead(track);
        }
        if remove {
            self.movie.movie.remove_track(track);
            self.movie.selected = None;
            self.movie.touch();
        }
    }

    fn track_lane(&mut self, ui: &mut Ui, track: usize, row: Rect, axis: TimeAxis) {
        let t = Tokens::current(ui.ctx());
        let lane = Rect::from_min_max(
            Pos2::new(axis.left, row.top() + 2.0),
            Pos2::new(axis.x(self.movie.movie.duration), row.bottom() - 2.0),
        );
        ui.painter()
            .rect_filled(lane, CornerRadius::same(radius::CONTROL), t.bg);
        let clips = self.movie.movie.tracks[track].clips.clone();
        for (index, clip) in clips.iter().enumerate() {
            self.clip_bar(ui, (track, index), clip, lane, axis);
        }
    }

    fn clip_bar(
        &mut self,
        ui: &mut Ui,
        at: (usize, usize),
        clip: &Clip,
        lane: Rect,
        axis: TimeAxis,
    ) {
        let (track, index) = at;
        let t = Tokens::current(ui.ctx());
        let bar = Rect::from_min_max(
            Pos2::new(axis.x(clip.start), lane.top()),
            Pos2::new(axis.x(clip.end()), lane.bottom()),
        );
        let selected = self.movie.selected == Some(at);
        let kind = self.movie.movie.tracks[track].kind;
        paint_bar(ui, bar, clip, clip_color(kind, &t), selected, &t);
        for (grab, area) in grab_areas(bar) {
            let response = ui.interact(
                area,
                Id::new(("movie-clip", track, index, grab as u8)),
                Sense::click_and_drag(),
            );
            self.handle_clip_response(&response, at, grab, axis);
        }
    }

    fn handle_clip_response(
        &mut self,
        response: &Response,
        at: (usize, usize),
        grab: Grab,
        axis: TimeAxis,
    ) {
        let (track, index) = at;
        if response.hovered() || response.dragged() {
            response_cursor(response, grab);
        }
        if response.clicked() || response.drag_started() {
            self.movie.selected = Some(at);
        }
        if response.dragged() {
            let dt = response.drag_delta().x / axis.pixels_per_second;
            if let Err(e) = self.movie.movie.drag_clip(track, index, grab, dt) {
                self.fail(e);
            }
            self.movie.touch();
        }
    }
}

fn response_cursor(response: &Response, grab: Grab) {
    let icon = match grab {
        Grab::Body => egui::CursorIcon::Grab,
        Grab::Left | Grab::Right => egui::CursorIcon::ResizeHorizontal,
    };
    response.ctx.set_cursor_icon(icon);
}

/// The resize handles at each end, then the body between them (later
/// interactions win, so the handles come last).
fn grab_areas(bar: Rect) -> [(Grab, Rect); 3] {
    let handle = HANDLE_WIDTH.min(bar.width() / 3.0);
    [
        (Grab::Body, bar),
        (
            Grab::Left,
            Rect::from_min_size(bar.min, Vec2::new(handle, bar.height())),
        ),
        (
            Grab::Right,
            Rect::from_min_size(
                Pos2::new(bar.right() - handle, bar.top()),
                Vec2::new(handle, bar.height()),
            ),
        ),
    ]
}

fn paint_bar(ui: &Ui, bar: Rect, clip: &Clip, color: Color32, selected: bool, t: &Tokens) {
    let painter = ui.painter();
    let corner = CornerRadius::same(radius::CONTROL);
    painter.rect_filled(bar, corner, color.gamma_multiply(0.85));
    if selected {
        painter.rect_stroke(bar, corner, Stroke::new(2.0, t.text), StrokeKind::Inside);
    }
    painter.with_clip_rect(bar).text(
        bar.left_center() + Vec2::new(space::TIGHT + 2.0, 0.0),
        Align2::LEFT_CENTER,
        clip.effect.label(),
        FontId::proportional(text::CAPTION),
        t.on_primary,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_axis_maps_seconds_to_x_and_back() {
        let axis = TimeAxis::new(
            Rect::from_min_size(Pos2::new(100.0, 0.0), Vec2::new(500.0, 10.0)),
            10.0,
        );
        assert_eq!(axis.x(0.0), 100.0);
        assert_eq!(axis.x(10.0), 600.0);
        assert_eq!(axis.seconds(350.0), 5.0);
    }

    #[test]
    fn ticks_thin_out_as_the_movie_gets_longer() {
        assert_eq!(tick_step(500.0), 0.25);
        assert_eq!(tick_step(10.0), 10.0);
        assert_eq!(tick_step(0.001), 300.0);
    }

    #[test]
    fn the_end_handles_never_swallow_a_short_bar() {
        let bar = Rect::from_min_size(Pos2::ZERO, Vec2::new(9.0, 20.0));
        let [_, left, right] = grab_areas(bar);
        assert!(left.1.right() <= right.1.left());
    }
}
