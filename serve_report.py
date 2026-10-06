"""PDAC 리포트 서버 - 생성된 HTML 리포트를 이 서버에서 그대로 띄운다.

파이프라인(run_pdac_agent.py)과 완전히 분리되어 있고, 표준 라이브러리만 씁니다.
biomni 패키지도, 가상환경도, 외부 서비스 로그인도 필요하지 않습니다.

    python serve_report.py                      # ./pdac_run 을 127.0.0.1:8000 에 서비스
    python serve_report.py --port 8080          # 포트 변경
    python serve_report.py --host 0.0.0.0       # 같은 네트워크의 다른 PC에서도 접근
    python serve_report.py --background         # 백그라운드로 띄우고 바로 프롬프트 복귀
    python serve_report.py --status             # 실행 여부 확인
    python serve_report.py --stop               # 종료

`/` 로 접속하면 디렉토리에 있는 리포트와 그림 목록이 나옵니다.
"""

import argparse
import http.server
import os
import signal
import socket
import socketserver
import subprocess
import sys
import time

DEFAULT_DIR = "./pdac_run"
DEFAULT_PORT = 8000
PID_FILE = ".report_server.pid"
LOG_FILE = "report_server.log"

REPORT_NAMES = ("pdac_report.html", "report.html", "index.html")


def _pid_path(directory: str) -> str:
    return os.path.join(directory, PID_FILE)


def _running_pid(directory: str) -> int | None:
    """Return the PID of a server already serving this directory, or None."""
    path = _pid_path(directory)
    if not os.path.exists(path):
        return None
    try:
        pid = int(open(path).read().strip())
        os.kill(pid, 0)  # signal 0 = existence check only
        return pid
    except (ValueError, ProcessLookupError, PermissionError, OSError):
        # Stale PID file from a crashed or killed server.
        os.remove(path)
        return None


def _local_ip() -> str | None:
    try:
        return socket.gethostbyname(socket.gethostname())
    except OSError:
        return None


def _entry_point(directory: str) -> str | None:
    for name in REPORT_NAMES:
        if os.path.exists(os.path.join(directory, name)):
            return name
    return None


def _index_html(directory: str) -> bytes:
    """A plain listing of the reports and figures in the directory."""
    reports, figures, data = [], [], []
    for name in sorted(os.listdir(directory)):
        if name.startswith("."):
            continue
        lower = name.lower()
        if lower.endswith((".html", ".htm")):
            reports.append(name)
        elif lower.endswith((".png", ".jpg", ".jpeg", ".svg")):
            figures.append(name)
        elif lower.endswith((".csv", ".tsv", ".json")):
            data.append(name)

    def section(title, names):
        if not names:
            return ""
        items = "".join(
            f'<li><a href="{name}">{name}</a> '
            f'<span class="size">{os.path.getsize(os.path.join(directory, name)) / 1024:.0f} KB</span></li>'
            for name in names
        )
        return f"<h2>{title}</h2><ul>{items}</ul>"

    body = (
        "<!doctype html><html lang='ko'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width, initial-scale=1'>"
        "<title>PDAC 리포트</title><style>"
        "body{font:15px/1.6 system-ui,-apple-system,'Segoe UI',sans-serif;max-width:720px;"
        "margin:48px auto;padding:0 24px;color:#0e1211;background:#f4f6f5}"
        "h1{font-size:22px;margin:0 0 4px}h2{font-size:13px;text-transform:uppercase;letter-spacing:.08em;"
        "color:#828b88;margin:28px 0 8px}ul{list-style:none;padding:0;margin:0}"
        "li{padding:8px 0;border-bottom:1px solid #dde2e0;display:flex;justify-content:space-between}"
        "a{color:#1c5cab;text-decoration:none}a:hover{text-decoration:underline}"
        ".size{font-family:ui-monospace,monospace;font-size:12px;color:#828b88}"
        ".path{font-family:ui-monospace,monospace;font-size:12px;color:#828b88;margin-bottom:8px}"
        "@media(prefers-color-scheme:dark){body{background:#0e100f;color:#fff}"
        "li{border-color:#2c2e2c}a{color:#6da7ec}}"
        "</style></head><body>"
        f"<h1>PDAC 리포트</h1><div class='path'>{os.path.abspath(directory)}</div>"
        + section("리포트", reports)
        + section("그림", figures)
        + section("데이터", data)
        + "</body></html>"
    )
    return body.encode("utf-8")


def make_handler(directory: str):
    class Handler(http.server.SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=directory, **kwargs)

        def do_GET(self):
            if self.path in ("/", "/index.html") and not os.path.exists(os.path.join(directory, "index.html")):
                payload = _index_html(directory)
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.send_header("Content-Length", str(len(payload)))
                self.end_headers()
                self.wfile.write(payload)
                return
            super().do_GET()

        def log_message(self, *_args):
            pass  # keep the console readable; the access log adds nothing here

    return Handler


def start_foreground(directory: str, host: str, port: int) -> None:
    socketserver.TCPServer.allow_reuse_address = True
    try:
        server = socketserver.TCPServer((host, port), make_handler(directory))
    except OSError as error:
        raise SystemExit(f"포트 {port} 를 열 수 없습니다 ({error}). 다른 포트를 쓰세요: --port 8080") from error

    entry = _entry_point(directory)
    print("=" * 70)
    print(f"리포트 서버 - {os.path.abspath(directory)}")
    print("=" * 70)
    print(f"  http://localhost:{port}/" + (entry or ""))
    if host == "0.0.0.0":
        ip = _local_ip()
        if ip:
            print(f"  http://{ip}:{port}/{entry or ''}   (같은 네트워크의 다른 PC)")
        print("  주의: 0.0.0.0 은 이 서버에 접근 가능한 모두에게 리포트를 공개합니다.")
    else:
        print("  (이 서버 안에서만 접근 가능. 다른 PC에서 열려면 --host 0.0.0.0)")
    if not entry:
        print(f"\n  경고: {directory} 에 HTML 리포트가 없습니다. 먼저 생성하세요:")
        print("        python run_pdac_agent.py --mode direct --stages 1-6")
    print("\nVS Code 원격 세션이면 포트가 자동 전달됩니다 - PORTS 탭의 지구본 아이콘으로 열 수 있습니다.")
    print("Ctrl+C 로 종료")

    with open(_pid_path(directory), "w") as handle:
        handle.write(str(os.getpid()))
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\n서버를 종료했습니다.")
    finally:
        server.server_close()
        if os.path.exists(_pid_path(directory)):
            os.remove(_pid_path(directory))


def start_background(directory: str, host: str, port: int) -> None:
    existing = _running_pid(directory)
    if existing:
        raise SystemExit(f"이미 PID {existing} 로 실행 중입니다. 종료하려면: python serve_report.py --stop")

    log_path = os.path.join(directory, LOG_FILE)
    with open(log_path, "a") as log:
        process = subprocess.Popen(
            [sys.executable, os.path.abspath(__file__), "--dir", directory, "--host", host, "--port", str(port)],
            stdout=log, stderr=log, start_new_session=True,
        )

    # Verify the child is ALIVE, not merely that something answers on the port: when another service
    # already holds the port, a connection probe succeeds while our own server is dead on arrival.
    deadline = time.time() + 5
    bound = False
    while time.time() < deadline:
        if process.poll() is not None:
            bound = False
            break
        with socket.socket() as probe:
            probe.settimeout(0.3)
            if probe.connect_ex(("127.0.0.1" if host == "0.0.0.0" else host, port)) == 0:
                time.sleep(0.4)  # a bind failure exits within a few hundred ms of the probe
                bound = process.poll() is None
                break
        time.sleep(0.2)

    if not bound:
        process.terminate()
        tail = ""
        if os.path.exists(log_path):
            with open(log_path) as log:
                tail = "".join(log.readlines()[-3:]).strip()
        raise SystemExit(f"서버가 시작되지 않았습니다.\n{tail}")

    with open(_pid_path(directory), "w") as handle:
        handle.write(str(process.pid))
    entry = _entry_point(directory) or ""
    print(f"백그라운드로 실행했습니다 (PID {process.pid})")
    print(f"  http://localhost:{port}/{entry}")
    print(f"  로그: {log_path}")
    print("  종료: python serve_report.py --stop")


def stop(directory: str) -> None:
    pid = _running_pid(directory)
    if pid is None:
        print("실행 중인 리포트 서버가 없습니다.")
        return
    os.kill(pid, signal.SIGTERM)
    if os.path.exists(_pid_path(directory)):
        os.remove(_pid_path(directory))
    print(f"PID {pid} 서버를 종료했습니다.")


def status(directory: str, port: int) -> None:
    pid = _running_pid(directory)
    if pid is None:
        print(f"실행 중이 아닙니다 ({os.path.abspath(directory)})")
        return
    entry = _entry_point(directory) or ""
    print(f"실행 중 - PID {pid}")
    print(f"  http://localhost:{port}/{entry}")


def main() -> None:
    parser = argparse.ArgumentParser(description="PDAC 리포트를 이 서버에서 HTTP 로 서비스한다")
    parser.add_argument("--dir", default=DEFAULT_DIR, help=f"리포트가 있는 디렉토리 (기본 {DEFAULT_DIR})")
    parser.add_argument("--port", type=int, default=DEFAULT_PORT, help=f"포트 (기본 {DEFAULT_PORT})")
    parser.add_argument("--host", default="127.0.0.1", help="바인딩 주소. 다른 PC에서 접근하려면 0.0.0.0")
    parser.add_argument("--background", action="store_true", help="백그라운드로 실행하고 프롬프트로 복귀")
    parser.add_argument("--stop", action="store_true", help="백그라운드 서버 종료")
    parser.add_argument("--status", action="store_true", help="실행 여부 확인")
    args = parser.parse_args()

    directory = os.path.abspath(args.dir)
    if not os.path.isdir(directory):
        raise SystemExit(f"{directory} 디렉토리가 없습니다. --dir 로 리포트 위치를 지정하세요.")

    if args.stop:
        stop(directory)
    elif args.status:
        status(directory, args.port)
    elif args.background:
        start_background(directory, args.host, args.port)
    else:
        start_foreground(directory, args.host, args.port)


if __name__ == "__main__":
    main()
